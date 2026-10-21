/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#include "eon/ChemfilesIO.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/ConFileIO.h"
#include "eon/HelperFunctions.h"
#include "eon/Matter.h"
#include "eon/Parameters.h"
#include "eon/ParametersJSON.h"

#include <chrono>
#include <filesystem>
#include <format>
#include <fstream>
#include <sstream>
#include <string>

#ifndef _WIN32
#include <fcntl.h>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
extern char **environ;
#else
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#endif

namespace {

namespace fs = std::filesystem;
using eonc::io::IoStatus;

static eonc::helpers::test::QuillTestLogger _quill_setup;

fs::path scratch(std::string_view tag) {
  const auto stamp =
      std::chrono::high_resolution_clock::now().time_since_epoch().count();
  auto dir = fs::temp_directory_path() / std::format("chfl_{}_{}", tag, stamp);
  fs::create_directories(dir);
  return dir;
}

void write_text(const fs::path &p, const std::string &text) {
  std::ofstream out(p);
  out << text;
}

std::string read_all(const fs::path &path) {
  std::ifstream in(path);
  std::stringstream buffer;
  buffer << in.rdbuf();
  return buffer.str();
}

// Spawn eonclient in `dir` with the test environment intact. std::system goes
// through a shell: macOS strips DYLD_* across /bin/sh, and cmd.exe quoting
// on Windows can run the client in the wrong directory (exit 1, log lost).
int run_eonclient(const fs::path &dir, const char *exe) {
  const fs::path log_path = dir / "client.log";
#ifdef _WIN32
  SECURITY_ATTRIBUTES sa{};
  sa.nLength = sizeof(sa);
  sa.bInheritHandle = TRUE;
  HANDLE log = CreateFileW(log_path.c_str(), GENERIC_WRITE, FILE_SHARE_READ,
                           &sa, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
  if (log == INVALID_HANDLE_VALUE) {
    return -1;
  }
  STARTUPINFOW si{};
  si.cb = sizeof(si);
  si.dwFlags = STARTF_USESTDHANDLES;
  si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
  si.hStdOutput = log;
  si.hStdError = log;
  std::wstring cmdline = L"\"" + fs::path(exe).wstring() + L"\"";
  PROCESS_INFORMATION pi{};
  const BOOL started = CreateProcessW(nullptr, cmdline.data(), nullptr, nullptr,
                                      TRUE, 0, nullptr, dir.c_str(), &si, &pi);
  CloseHandle(log);
  if (!started) {
    return -1;
  }
  WaitForSingleObject(pi.hProcess, INFINITE);
  DWORD code = 1;
  GetExitCodeProcess(pi.hProcess, &code);
  CloseHandle(pi.hThread);
  CloseHandle(pi.hProcess);
  return static_cast<int>(code);
#else
  const std::string log_name = log_path.string();
  posix_spawn_file_actions_t actions;
  posix_spawn_file_actions_init(&actions);
  posix_spawn_file_actions_addopen(&actions, STDOUT_FILENO, log_name.c_str(),
                                   O_WRONLY | O_CREAT | O_TRUNC, 0644);
  posix_spawn_file_actions_adddup2(&actions, STDOUT_FILENO, STDERR_FILENO);
  char *argv[] = {const_cast<char *>(exe), nullptr};
  pid_t pid = 0;
  int spawned = -1;
  int status = 1;
  {
    const fs::path previous = fs::current_path();
    fs::current_path(dir);
    spawned = posix_spawn(&pid, exe, &actions, nullptr, argv, environ);
    std::error_code ec;
    fs::current_path(previous, ec);
  }
  posix_spawn_file_actions_destroy(&actions);
  if (spawned != 0) {
    return -1;
  }
  if (waitpid(pid, &status, 0) < 0) {
    return -1;
  }
  if (WIFEXITED(status)) {
    return WEXITSTATUS(status);
  }
  if (WIFSIGNALED(status)) {
    return 128 + WTERMSIG(status);
  }
  return 1;
#endif
}

/// Three Ar atoms near the LJ minimum spacing, with a 15 A cubic lattice.
const char *kArLattice = "3\n"
                         "Lattice=\"15.0 0.0 0.0 0.0 15.0 0.0 0.0 0.0 15.0\" "
                         "Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n"
                         "Ar 5.0 5.0 5.0\n"
                         "Ar 6.2 5.0 5.0\n"
                         "Ar 5.6 6.0 5.0\n";

/// Same atoms, no cell at all.
const char *kArNoCell = "3\n"
                        "plain xyz comment\n"
                        "Ar 1.0 2.0 3.0\n"
                        "Ar 2.2 2.0 3.0\n"
                        "Ar 1.6 3.0 3.0\n";

struct Fixture {
  Parameters params;
  std::shared_ptr<Potential> pot;
  Fixture() {
    ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
    pot = eonc::helpers::sharePotential(
        eonc::helpers::makePotential(PotType::LJ, params));
    eonc::io::set_config_cell({}, true);
    eonc::io::set_default_vacuum(10.0);
    eonc::io::set_structure_input_format("");
  }
  ~Fixture() {
    eonc::io::set_config_cell({}, true);
    eonc::io::set_default_vacuum(10.0);
    eonc::io::set_structure_input_format("");
  }
  std::shared_ptr<Matter> matter() {
    return std::make_shared<Matter>(pot, params);
  }
};

} // namespace

TEST_CASE("Path routing picks readcon or chemfiles", "[chemfiles][io]") {
  using eonc::io::is_chemfiles_path;
  using eonc::io::is_readcon_path;
  CHECK(is_readcon_path("pos.con"));
  CHECK(is_readcon_path("band.con.gz"));
  CHECK(is_readcon_path("pos.convel"));
  CHECK_FALSE(is_readcon_path("pos.xyz"));
  CHECK(is_chemfiles_path("pos.xyz"));
  CHECK(is_chemfiles_path("pos.extxyz"));
  CHECK(is_chemfiles_path("x.PDB"));
  CHECK(is_chemfiles_path("POSCAR"));
  CHECK_FALSE(is_chemfiles_path("pos.con"));
  CHECK_FALSE(is_chemfiles_path("pos_0001"));
}

TEST_CASE("[Main] cell parsing", "[chemfiles][config]") {
  using eonc::io::parse_cell_values;
  CHECK(parse_cell_values("").empty());
  CHECK(parse_cell_values("10 11 12").size() == 3);
  CHECK(parse_cell_values("10 0 0  0 10 0  0 0 10").size() == 9);
  CHECK_THROWS_AS(parse_cell_values("10 11"), std::invalid_argument);
  CHECK_THROWS_AS(parse_cell_values("10 x 12"), std::invalid_argument);
}

TEST_CASE_METHOD(Fixture, "inp2matter still reads .con through readcon",
                 "[chemfiles][io]") {
  auto m = matter();
  REQUIRE(eonc::io::io_ok(m->inp2matter("reactant.con")));
  REQUIRE(m->numberOfAtoms() > 0);
}

TEST_CASE_METHOD(Fixture, "Foreign file without chemfiles is refused clearly",
                 "[chemfiles][io]") {
  if (eonc::io::chemfiles_available()) {
    SKIP("chemfiles linked");
  }
  const auto dir = scratch("nochfl");
  write_text(dir / "pos.xyz", kArLattice);
  auto m = matter();
  CHECK(m->inp2matter((dir / "pos.xyz").string()) == IoStatus::InvalidArgument);
}

TEST_CASE_METHOD(Fixture, "Cell path 1: the file's own lattice wins",
                 "[chemfiles][cell]") {
  if (!eonc::io::chemfiles_available()) {
    SKIP("built without chemfiles");
  }
  // A config cell must not override a cell the file carries.
  eonc::io::set_config_cell({30.0, 30.0, 30.0}, false);
  const auto dir = scratch("cell1");
  write_text(dir / "pos.xyz", kArLattice);
  auto m = matter();
  REQUIRE(eonc::io::io_ok(m->inp2matter((dir / "pos.xyz").string())));
  REQUIRE(m->numberOfAtoms() == 3);
  const Matrix3d cell = m->getCell();
  CHECK(cell(0, 0) == Catch::Approx(15.0));
  CHECK(cell(1, 1) == Catch::Approx(15.0));
  CHECK(cell(2, 2) == Catch::Approx(15.0));
  CHECK(cell(0, 1) == Catch::Approx(0.0).margin(1e-12));
  CHECK(m->getPeriodic());
  CHECK(m->getAtomicNr(0) == 18);
  CHECK(m->getMass(0) == Catch::Approx(39.948).margin(0.01));
  CHECK(m->getPosition(1, 0) == Catch::Approx(6.2));
  CHECK(eonc::io::last_cell_warning().empty());
}

TEST_CASE_METHOD(Fixture,
                 "Cell path 2: [Main] cell and pbc fill a missing cell",
                 "[chemfiles][cell]") {
  if (!eonc::io::chemfiles_available()) {
    SKIP("built without chemfiles");
  }
  const auto dir = scratch("cell2");
  write_text(dir / "pos.xyz", kArNoCell);

  SECTION("three lengths") {
    eonc::io::set_config_cell({20.0, 21.0, 22.0}, true);
    auto m = matter();
    REQUIRE(eonc::io::io_ok(m->inp2matter((dir / "pos.xyz").string())));
    const Matrix3d cell = m->getCell();
    CHECK(cell(0, 0) == Catch::Approx(20.0));
    CHECK(cell(1, 1) == Catch::Approx(21.0));
    CHECK(cell(2, 2) == Catch::Approx(22.0));
    CHECK(m->getPeriodic());
    CHECK(m->getPosition(0, 0) == Catch::Approx(1.0)); // not shifted
    CHECK(eonc::io::last_cell_warning().empty());
  }
  SECTION("full lattice, pbc off") {
    eonc::io::set_config_cell({20, 0, 0, 2, 20, 0, 0, 0, 25}, false);
    auto m = matter();
    REQUIRE(eonc::io::io_ok(m->inp2matter((dir / "pos.xyz").string())));
    const Matrix3d cell = m->getCell();
    CHECK(cell(1, 0) == Catch::Approx(2.0));
    CHECK(cell(2, 2) == Catch::Approx(25.0));
    CHECK_FALSE(m->getPeriodic());
  }
}

TEST_CASE_METHOD(Fixture,
                 "Cell path 3: no cell anywhere warns and pads the atoms",
                 "[chemfiles][cell]") {
  if (!eonc::io::chemfiles_available()) {
    SKIP("built without chemfiles");
  }
  const auto dir = scratch("cell3");
  const auto file = (dir / "pos.xyz").string();
  write_text(dir / "pos.xyz", kArNoCell);
  eonc::io::set_default_vacuum(7.5);
  auto m = matter();
  REQUIRE(eonc::io::io_ok(m->inp2matter(file)));

  const Matrix3d cell = m->getCell();
  CHECK(cell(0, 0) == Catch::Approx(1.2 + 15.0));
  CHECK(cell(1, 1) == Catch::Approx(1.0 + 15.0));
  CHECK(cell(2, 2) == Catch::Approx(0.0 + 15.0));
  CHECK(cell(0, 1) == Catch::Approx(0.0).margin(1e-12));
  CHECK_FALSE(m->getPeriodic());
  // Lowest corner of the atoms sits vacuum away from the origin.
  CHECK(m->getPosition(0, 0) == Catch::Approx(7.5));
  CHECK(m->getPosition(0, 1) == Catch::Approx(7.5));
  CHECK(m->getPosition(0, 2) == Catch::Approx(7.5));

  const auto &warn = eonc::io::last_cell_warning();
  CHECK(warn.starts_with("WARNING"));
  CHECK(warn.find(file) != std::string::npos);
  CHECK(warn.find("[Main] cell") != std::string::npos);
  CHECK(warn.find("[Main] default_vacuum") != std::string::npos);
}

TEST_CASE_METHOD(Fixture, "extxyz fixed atoms and velocities map into Matter",
                 "[chemfiles][fixed]") {
  if (!eonc::io::chemfiles_available()) {
    SKIP("built without chemfiles");
  }
  const auto dir = scratch("fixed");
  SECTION("ASE move_mask, true means free") {
    write_text(dir / "mm.extxyz",
               "3\n"
               "Lattice=\"15 0 0 0 15 0 0 0 15\" "
               "Properties=species:S:1:pos:R:3:move_mask:L:1\n"
               "Ar 5.0 5.0 5.0 F\n"
               "Ar 6.2 5.0 5.0 T\n"
               "Ar 5.6 6.0 5.0 T\n");
    auto m = matter();
    REQUIRE(eonc::io::io_ok(m->inp2matter((dir / "mm.extxyz").string())));
    CHECK(m->getFixed(0) != 0);
    CHECK(m->getFixed(1) == 0);
    CHECK(m->getFixed(2) == 0);
    CHECK(m->numberOfFixedAtoms() == 1);
  }
  SECTION("fixed as a per-axis vector") {
    write_text(dir / "fx.extxyz", "2\n"
                                  "Lattice=\"15 0 0 0 15 0 0 0 15\" "
                                  "Properties=species:S:1:pos:R:3:fixed:R:3\n"
                                  "Ar 5.0 5.0 5.0 0 0 1\n"
                                  "Ar 6.2 5.0 5.0 0 0 0\n");
    auto m = matter();
    REQUIRE(eonc::io::io_ok(m->inp2matter((dir / "fx.extxyz").string())));
    CHECK(m->getFixed(0, 0) == 0);
    CHECK(m->getFixed(0, 2) != 0);
    CHECK(m->getFixed(1, 2) == 0);
  }
  SECTION("velocities") {
    write_text(dir / "v.extxyz", "2\n"
                                 "Lattice=\"15 0 0 0 15 0 0 0 15\" "
                                 "Properties=species:S:1:pos:R:3:velo:R:3\n"
                                 "Ar 5.0 5.0 5.0 0.1 0.2 0.3\n"
                                 "Ar 6.2 5.0 5.0 0.0 0.0 -0.1\n");
    auto m = matter();
    REQUIRE(eonc::io::io_ok(m->inp2matter((dir / "v.extxyz").string())));
    const AtomMatrix v = m->getVelocities();
    CHECK(v(0, 1) == Catch::Approx(0.2));
    CHECK(v(1, 2) == Catch::Approx(-0.1));
  }
}

TEST_CASE_METHOD(Fixture, "xyz, extxyz and pdb round-trip through .con",
                 "[chemfiles][roundtrip]") {
  if (!eonc::io::chemfiles_available()) {
    SKIP("built without chemfiles");
  }
  const auto dir = scratch("rt");
  write_text(dir / "a.xyz", kArLattice);
  write_text(dir / "b.extxyz", "3\n"
                               "Lattice=\"15 0 0 0 15 0 0 0 15\" "
                               "Properties=species:S:1:pos:R:3:move_mask:L:1\n"
                               "Ar 5.0 5.0 5.0 F\n"
                               "Ar 6.2 5.0 5.0 T\n"
                               "Ar 5.6 6.0 5.0 T\n");
  write_text(dir / "c.pdb",
             "CRYST1   15.000   15.000   15.000  90.00  90.00  90.00 P 1   "
             "        1\n"
             "HETATM    1 AR    AR A   1       5.000   5.000   5.000  1.00  "
             "0.00          AR\n"
             "HETATM    2 AR    AR A   1       6.200   5.000   5.000  1.00  "
             "0.00          AR\n"
             "HETATM    3 AR    AR A   1       5.600   6.000   5.000  1.00  "
             "0.00          AR\n"
             "END\n");
  for (const char *name : {"a.xyz", "b.extxyz", "c.pdb"}) {
    INFO(name);
    auto in = matter();
    REQUIRE(eonc::io::io_ok(in->inp2matter((dir / name).string())));
    const auto con = (dir / (std::string(name) + ".con")).string();
    REQUIRE(eonc::io::io_ok(eonc::io::matter2con(*in, con)));
    auto back = matter();
    REQUIRE(eonc::io::io_ok(back->inp2matter(con)));
    REQUIRE(back->numberOfAtoms() == in->numberOfAtoms());
    CHECK((back->getPositions() - in->getPositions()).cwiseAbs().maxCoeff() <
          1e-3);
    CHECK((back->getCell() - in->getCell()).cwiseAbs().maxCoeff() < 1e-6);
    CHECK(back->getAtomicNrs() == in->getAtomicNrs());
    CHECK(back->numberOfFixedAtoms() == in->numberOfFixedAtoms());
  }
}

TEST_CASE_METHOD(Fixture, "input_format lets pos.con resolve to pos.xyz",
                 "[chemfiles][io]") {
  if (!eonc::io::chemfiles_available()) {
    SKIP("built without chemfiles");
  }
  const auto dir = scratch("fmt");
  write_text(dir / "pos.xyz", kArLattice);
  eonc::io::set_structure_input_format("xyz");
  auto m = matter();
  REQUIRE(eonc::io::io_ok(m->con2matter((dir / "pos.con").string())));
  CHECK(m->numberOfAtoms() == 3);
}

TEST_CASE("eonclient minimizes and evaluates a point from pos.xyz",
          "[chemfiles][e2e]") {
  if (!eonc::io::chemfiles_available()) {
    SKIP("built without chemfiles");
  }
  const char *exe = std::getenv("EONCLIENT");
  if (exe == nullptr) {
    SKIP("EONCLIENT not set");
  }
  for (const char *job : {"minimization", "point"}) {
    INFO(job);
    const auto dir = scratch(job);
    write_text(dir / "pos.xyz", kArLattice);
    write_text(dir / "config.ini", std::format("[Main]\n"
                                               "job = {}\n"
                                               "input_format = xyz\n"
                                               "\n"
                                               "[Potential]\n"
                                               "potential = lj\n",
                                               job));
    const int rc = run_eonclient(dir, exe);
    if (rc != 0) {
      INFO(read_all(dir / "client.log"));
    }
    REQUIRE(rc == 0);
    if (std::string(job) == "minimization") {
      CHECK(fs::exists(dir / "min.con"));
    } else {
      CHECK(fs::exists(dir / "results.dat"));
    }
  }
}

TEST_CASE("Main structure input is stored and published",
          "[chemfiles][config]") {
  const std::string format = eonc::io::structure_input_format();
  const double vacuum = eonc::io::default_vacuum();
  struct Restore {
    std::string format;
    double vacuum;
    ~Restore() {
      eonc::io::set_structure_input_format(format);
      eonc::io::set_config_cell({}, true);
      eonc::io::set_default_vacuum(vacuum);
    }
  } restore{format, vacuum};

  Parameters loaded;
  REQUIRE(loaded.load_ini_text("[Main]\n"
                               "input_format = xyz\n"
                               "cell = 20 21 22\n"
                               "pbc = false\n"
                               "default_vacuum = 7.5\n") == 0);
  REQUIRE(loaded.main_options().inputFormat == "xyz");
  REQUIRE(loaded.main_options().cell == "20 21 22");
  REQUIRE_FALSE(loaded.main_options().pbc);
  REQUIRE(loaded.main_options().defaultVacuum == Catch::Approx(7.5));
  REQUIRE(eonc::io::structure_input_format() == "xyz");
  REQUIRE(eonc::io::default_vacuum() == Catch::Approx(7.5));

  eonc::io::set_structure_input_format("");
  eonc::io::set_default_vacuum(10.0);
  Parameters round;
  eonc::config::from_json(eonc::config::to_json(loaded), round);
  REQUIRE(round.main_options().inputFormat == "xyz");
  REQUIRE(round.main_options().cell == "20 21 22");
  REQUIRE_FALSE(round.main_options().pbc);
  REQUIRE(round.main_options().defaultVacuum == Catch::Approx(7.5));
  REQUIRE(eonc::io::structure_input_format() == "xyz");
  REQUIRE(eonc::io::default_vacuum() == Catch::Approx(7.5));

  Parameters con;
  REQUIRE(con.load_ini_text("[Main]\ninput_format = con\n") == 0);
  REQUIRE(con.main_options().inputFormat == "con");
  REQUIRE(eonc::io::structure_input_format().empty());

  Parameters bad;
  REQUIRE_THROWS_AS(bad.load_ini_text("[Main]\ncell = 1 2\n"),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(bad.load_ini_text("[Main]\ndefault_vacuum = -1\n"),
                    std::invalid_argument);
  REQUIRE(bad.main_options().cell.empty());
  REQUIRE(bad.main_options().defaultVacuum == Catch::Approx(10.0));
}
