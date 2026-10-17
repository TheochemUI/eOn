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
#include "ReadconDbMirror.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/ConFileIO.h"
#include "eon/HelperFunctions.h"
#include "eon/Matter.h"
#include "eon/Parameters.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace {

namespace fs = std::filesystem;

static eonc::helpers::test::QuillTestLogger _quill_setup;

std::uint64_t bits_of(double value) {
  std::uint64_t out = 0;
  std::memcpy(&out, &value, sizeof(out));
  return out;
}

double from_bits(std::uint64_t bits) {
  double out = 0.0;
  std::memcpy(&out, &bits, sizeof(out));
  return out;
}

std::string read_bytes(const fs::path &path) {
  std::ifstream in(path, std::ios::binary);
  REQUIRE(in.good());
  return std::string(std::istreambuf_iterator<char>(in),
                     std::istreambuf_iterator<char>());
}

fs::path golden_path(const char *name) {
  return fs::path("../../restart") / name;
}

Matrix3d cell_from_lengths_angles() {
  constexpr double kPi = 3.14159265358979323846;
  // Beta is one ulp above 75 degrees. Lengths 10, 11, 12 and these
  // angles are a fixed point of the cell map, so a read then write
  // keeps the same length and angle bytes.
  const double alpha = 80.0 * kPi / 180.0;
  const double beta = 0x1.2c00000000001p+6 * kPi / 180.0;
  const double gamma = 100.0 * kPi / 180.0;
  Matrix3d cell = Matrix3d::Zero();
  cell(0, 0) = 1.0;
  cell(1, 0) = std::cos(gamma);
  cell(1, 1) = std::sin(gamma);
  cell(2, 0) = std::cos(beta);
  cell(2, 1) = (std::cos(alpha) - cell(1, 0) * cell(2, 0)) / cell(1, 1);
  cell(2, 2) =
      std::sqrt(1.0 - cell(2, 0) * cell(2, 0) - cell(2, 1) * cell(2, 1));
  cell.row(0) *= 10.0;
  cell.row(1) *= 11.0;
  cell.row(2) *= 12.0;
  return cell;
}

struct MatterHolder {
  Parameters params;
  std::shared_ptr<Potential> pot;
  std::shared_ptr<Matter> matter;

  MatterHolder() {
    ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
    pot = eonc::helpers::sharePotential(
        eonc::helpers::makePotential(PotType::LJ, params));
    matter = std::make_shared<Matter>(pot, params);
  }
};

void fill_restart_matter(Matter &m, bool with_velocities) {
  m.resize(4);
  m.setPeriodic(false);
  m.setCell(cell_from_lengths_angles());
  const int z[4] = {29, 1, 29, 1};
  const double mass[4] = {63.546, 1.008, 63.546, 1.008};
  const double sub = from_bits(1);
  const double x[4] = {-0.0, 1.0, 2.0, 3.0};
  const double y[4] = {0.5, 1.5, 2.5, 3.5};
  const double zpos[4] = {sub, 1.5, 2.5, 3.5};
  for (long i = 0; i < 4; ++i) {
    m.setAtomicNr(i, z[i]);
    m.setMass(i, mass[i]);
    m.setAtomIndex(i, i + 1);
    m.setFixed(i, false);
    m.setPosition(i, 0, x[i]);
    m.setPosition(i, 1, y[i]);
    m.setPosition(i, 2, zpos[i]);
  }
  m.setFixed(0, 0, true);

  AtomMatrix forces(4, 3);
  forces << 0.0, 0.1, -0.2, 0.3, -0.0, 0.4, -0.5, 0.2, 0.1, 0.05, -0.3, 0.25;
  m.setEvaluation(forces, -1.25);

  if (with_velocities) {
    AtomMatrix vel(4, 3);
    vel << 0.0, -0.0, sub, 0.01, -0.02, 0.03, -0.04, 0.05, -0.06, 0.07, 0.08,
        -0.09;
    m.setVelocities(vel);
  }
}

struct CwdGuard {
  fs::path previous;
  explicit CwdGuard(const fs::path &dir)
      : previous(fs::current_path()) {
    fs::create_directories(dir);
    fs::current_path(dir);
  }
  ~CwdGuard() {
    std::error_code ec;
    fs::current_path(previous, ec);
  }
};

void write_run(Matter &m, int steps, bool resume) {
  for (int i = 0; i < steps; ++i) {
    auto pos = m.getPositions();
    pos(0, 0) += 0.01;
    m.setPositions(pos);
    REQUIRE(eonc::io::io_ok(m.matter2con("pos_cp.con", false)));
    const bool append = resume || i > 0 || fs::exists("movie.con");
    REQUIRE(eonc::io::io_ok(m.matter2con("movie.con", append)));
  }
}

} // namespace

TEST_CASE("corpus path id is FNV-1a of the path only", "[con][corpus]") {
  REQUIRE(eonc::io::corpus_path_id("/var/eon/restart_exact.con") ==
          14164577628502097180ull);
}

TEST_CASE("appended corpus frames flush in groups of 64", "[con][corpus]") {
#if defined(_WIN32)
  SUCCEED("the mirror is compiled out on Windows");
#else
  const auto dir =
      fs::temp_directory_path() / "eon-corpus-buffer-restart-exact";
  fs::remove_all(dir);
  fs::create_directories(dir);
  const auto path = (dir / "movie.con").string();
  eonc::io::flush_corpus_appends();
  const std::string frame = "frame\n";
  for (int i = 0; i < 63; ++i) {
    eonc::io::note_corpus_append(path, frame);
  }
  REQUIRE(eonc::io::corpus_buffered_frames(path) == 63);
  eonc::io::note_corpus_append(path, frame);
  REQUIRE(eonc::io::corpus_buffered_frames(path) == 0);
  eonc::io::note_corpus_append(path, frame);
  REQUIRE(eonc::io::corpus_buffered_frames(path) == 1);
  eonc::io::flush_corpus_appends();
  REQUIRE(eonc::io::corpus_buffered_frames(path) == 0);
  fs::remove_all(dir);
#endif
}

TEST_CASE("matter2con matches the restart golden", "[con][restart]") {
  MatterHolder hold;
  fill_restart_matter(*hold.matter, false);
  const auto tmp = fs::temp_directory_path() / "eon_restart_exact_out.con";
  eonc::io::ConFrameMetadata meta;
  meta.energy = -1.25;
  meta.write_con_forces = true;
  REQUIRE(eonc::io::io_ok(hold.matter->matter2con(tmp.string(), false, &meta)));
  const auto golden = golden_path("restart_exact.con");
  REQUIRE(read_bytes(tmp) == read_bytes(golden));

  MatterHolder back;
  back.matter->setPeriodic(false);
  REQUIRE(eonc::io::io_ok(back.matter->con2matter(golden.string())));
  const auto wrote = hold.matter->getPositions();
  const auto read = back.matter->getPositions();
  REQUIRE(wrote.rows() == read.rows());
  REQUIRE(std::memcmp(wrote.data(), read.data(),
                      static_cast<std::size_t>(wrote.size()) *
                          sizeof(double)) == 0);
  REQUIRE(bits_of(read(0, 0)) == bits_of(-0.0));
  REQUIRE(bits_of(read(0, 2)) == 1ull);

  eonc::io::set_write_con_forces(true);
  const auto again = fs::temp_directory_path() / "eon_restart_exact_again.con";
  REQUIRE(eonc::io::io_ok(back.matter->matter2con(again.string(), false)));
  eonc::io::set_write_con_forces(false);
  REQUIRE(read_bytes(again) == read_bytes(golden));
  fs::remove(tmp);
  fs::remove(again);
}

TEST_CASE("matter2convel matches the velocity golden", "[con][restart]") {
  MatterHolder hold;
  fill_restart_matter(*hold.matter, true);
  const auto tmp = fs::temp_directory_path() / "eon_restart_exact_vel.convel";
  REQUIRE(eonc::io::io_ok(eonc::io::matter2convel(*hold.matter, tmp.string())));
  const auto golden = golden_path("restart_exact_vel.con");
  REQUIRE(read_bytes(tmp) == read_bytes(golden));

  MatterHolder back;
  back.matter->setPeriodic(false);
  // The golden keeps a .con name. convel2matter would append another
  // .convel suffix. The frame is the convel writer output either way.
  REQUIRE(eonc::io::io_ok(back.matter->con2matter(golden.string())));
  const auto wrote = hold.matter->getVelocities();
  const auto read = back.matter->getVelocities();
  REQUIRE(std::memcmp(wrote.data(), read.data(),
                      static_cast<std::size_t>(wrote.size()) *
                          sizeof(double)) == 0);
  REQUIRE(bits_of(read(0, 1)) == bits_of(-0.0));
  REQUIRE(bits_of(read(0, 2)) == 1ull);
  fs::remove(tmp);
}

TEST_CASE("two Cu masses collapse to the first mass", "[con][mass]") {
  MatterHolder hold;
  auto &m = *hold.matter;
  m.resize(2);
  m.setPeriodic(false);
  m.setCell(Matrix3d::Identity() * 10.0);
  m.setAtomicNr(0, 29);
  m.setAtomicNr(1, 29);
  m.setMass(0, 63.5);
  m.setMass(1, 62.0);
  m.setAtomIndex(0, 1);
  m.setAtomIndex(1, 2);
  const auto tmp = fs::temp_directory_path() / "eon_two_cu_mass.con";
  REQUIRE(eonc::io::io_ok(m.matter2con(tmp.string(), false)));
  MatterHolder back;
  REQUIRE(eonc::io::io_ok(back.matter->con2matter(tmp.string())));
  REQUIRE(back.matter->getMass(0) == 63.5);
  REQUIRE(back.matter->getMass(1) == 63.5);
  fs::remove(tmp);
}

TEST_CASE("mid-run pos_cp restart matches an uninterrupted run",
          "[con][restart]") {
  const auto root = fs::temp_directory_path() / "eon-mid-run-restart";
  fs::remove_all(root);
  const auto full = root / "full";
  const auto stopped = root / "stopped";

  {
    MatterHolder hold;
    auto &m = *hold.matter;
    m.resize(2);
    m.setPeriodic(false);
    m.setCell(Matrix3d::Identity() * 20.0);
    m.setAtomicNr(0, 29);
    m.setAtomicNr(1, 1);
    m.setMass(0, 63.546);
    m.setMass(1, 1.008);
    m.setAtomIndex(0, 1);
    m.setAtomIndex(1, 2);
    CwdGuard guard(full);
    write_run(m, 4, false);
  }
  {
    MatterHolder hold;
    auto &m = *hold.matter;
    m.resize(2);
    m.setPeriodic(false);
    m.setCell(Matrix3d::Identity() * 20.0);
    m.setAtomicNr(0, 29);
    m.setAtomicNr(1, 1);
    m.setMass(0, 63.546);
    m.setMass(1, 1.008);
    m.setAtomIndex(0, 1);
    m.setAtomIndex(1, 2);
    CwdGuard guard(stopped);
    write_run(m, 2, false);
    eonc::io::resetConAppendState();
    MatterHolder resumed;
    auto &next = *resumed.matter;
    next.setPeriodic(false);
    REQUIRE(eonc::io::io_ok(next.con2matter("pos_cp.con")));
    write_run(next, 2, true);
  }

  REQUIRE(read_bytes(full / "movie.con") == read_bytes(stopped / "movie.con"));
  REQUIRE(read_bytes(full / "pos_cp.con") ==
          read_bytes(stopped / "pos_cp.con"));

  MatterHolder a;
  MatterHolder b;
  {
    CwdGuard guard(full);
    REQUIRE(eonc::io::io_ok(a.matter->con2matter("pos_cp.con")));
  }
  {
    CwdGuard guard(stopped);
    REQUIRE(eonc::io::io_ok(b.matter->con2matter("pos_cp.con")));
  }
  const auto pa = a.matter->getPositions();
  const auto pb = b.matter->getPositions();
  REQUIRE(std::memcmp(pa.data(), pb.data(),
                      static_cast<std::size_t>(pa.size()) * sizeof(double)) ==
          0);
  fs::remove_all(root);
}

TEST_CASE("saddle resume loads displacement_cp.con and mode_cp.dat",
          "[con][restart]") {
  const auto dir = fs::temp_directory_path() / "eon-saddle-resume";
  fs::remove_all(dir);
  MatterHolder hold;
  auto &m = *hold.matter;
  m.resize(2);
  m.setPeriodic(false);
  m.setCell(Matrix3d::Identity() * 15.0);
  m.setAtomicNr(0, 29);
  m.setAtomicNr(1, 1);
  m.setMass(0, 63.546);
  m.setMass(1, 1.008);
  m.setAtomIndex(0, 1);
  m.setAtomIndex(1, 2);
  m.setPosition(0, 0, 1.25);
  m.setPosition(1, 1, -0.5);
  AtomMatrix mode(2, 3);
  mode << 0.25, 0.0, -0.5, 0.0, 0.5, 0.25;

  {
    CwdGuard guard(dir);
    REQUIRE(eonc::io::io_ok(m.matter2con("displacement_cp.con", false)));
    eonc::helpers::saveMode("mode_cp.dat", hold.matter, mode);
  }

  MatterHolder back;
  back.matter->setPeriodic(false);
  AtomMatrix loaded(2, 3);
  {
    CwdGuard guard(dir);
    REQUIRE(eonc::io::io_ok(back.matter->con2matter("displacement_cp.con")));
    loaded = eonc::helpers::loadMode("mode_cp.dat", 2);
  }
  const auto wrote = m.getPositions();
  const auto read = back.matter->getPositions();
  REQUIRE(std::memcmp(wrote.data(), read.data(),
                      static_cast<std::size_t>(wrote.size()) *
                          sizeof(double)) == 0);
  REQUIRE(loaded(0, 0) == 0.25);
  REQUIRE(loaded(0, 2) == -0.5);
  REQUIRE(loaded(1, 1) == 0.5);
  REQUIRE(fs::exists(dir / "displacement_cp.con"));
  REQUIRE(fs::exists(dir / "mode_cp.dat"));
  fs::remove_all(dir);
}
