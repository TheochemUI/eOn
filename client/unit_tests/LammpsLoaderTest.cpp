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

/// Tests for the LAMMPS runtime loader.
/// These tests verify the dlopen-based loader interface works correctly
/// regardless of whether LAMMPS is actually installed.

#include "eon/potentials/LAMMPS/LammpsLoader.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include <cstdint>
#include <filesystem>
#include <fstream>
#include "eon/potentials/LAMMPS/LAMMPSPot.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#ifndef _WIN32
#include <unistd.h>
#endif

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

namespace {
class MockLammpsLoader : public eonc::ILammpsLoader {
public:
  int require_calls{0};
  void require_loaded() override { ++require_calls; }
  [[nodiscard]] bool is_loaded() const noexcept override { return true; }
  [[nodiscard]] bool available() const override { return true; }
  [[nodiscard]] const std::string &last_error() const noexcept override {
    return err_;
  }

private:
  std::string err_{};
};
} // namespace

#ifndef _WIN32
void *lammpsOpenStub(int, char **, void **) {
  return reinterpret_cast<void *>(static_cast<std::uintptr_t>(1));
}
void lammpsCloseStub(void *) {}
char *lammpsCommandStub(void *, const char *) { return nullptr; }
void lammpsScatterStub(void *, const char *, int, int, void *) {}
void *lammpsExtractNull(void *, const char *, const char *) { return nullptr; }

TEST_CASE("a missing LAMMPS variable rejects the geometry",
          "[lammps][worker]") {
  namespace fs = std::filesystem;
  const fs::path previous = fs::current_path();
  const fs::path dir = fs::temp_directory_path() / "eon-lmp-null";
  fs::create_directories(dir);
  fs::current_path(dir);
  {
    std::ofstream input("in.lammps");
    input << "pair_style none\n";
  }
  MockLammpsLoader mock;
  mock.open_no_mpi = lammpsOpenStub;
  mock.close = lammpsCloseStub;
  mock.command = lammpsCommandStub;
  mock.scatter_atoms = lammpsScatterStub;
  mock.extract_variable = lammpsExtractNull;
  eonc::Parameters params;
  LAMMPSPot pot(params, mock);
  const long n = 2;
  const double positions[6] = {0.0, 0.0, 0.0, 1.5, 0.0, 0.0};
  const int atomic[2] = {1, 1};
  double forces[6] = {};
  double energy = 0.0;
  const double box[9] = {10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0};
  pot.force(n, positions, atomic, forces, &energy, nullptr, box);
  REQUIRE(energy == Catch::Approx(1.0e6));
  REQUIRE(forces[0] == Catch::Approx(1.0));
  REQUIRE(forces[3] == Catch::Approx(-1.0));
  fs::current_path(previous);
  std::error_code ec;
  fs::remove_all(dir, ec);
}
#endif

TEST_CASE("LAMMPSPot constructs with injected loader mock",
          "[lammps][loader][inject]") {
  const bool singleton_loaded = eonc::LammpsLoader::instance().is_loaded();
  MockLammpsLoader mock;
  eonc::Parameters params;
  LAMMPSPot pot(params, mock);
  REQUIRE(mock.require_calls == 1);
  REQUIRE(eonc::LammpsLoader::instance().is_loaded() == singleton_loaded);
}

TEST_CASE("LammpsLoader: singleton returns consistent instance",
          "[lammps][loader]") {
  auto &a = eonc::LammpsLoader::instance();
  auto &b = eonc::LammpsLoader::instance();
  REQUIRE(&a == &b);
}

TEST_CASE("LammpsLoader: available probes without loading",
          "[lammps][loader]") {
  auto &loader = eonc::LammpsLoader::instance();
  // Filesystem probe must not dlopen, so it cannot move the loader from
  // unloaded to loaded. Probing twice must also agree with itself.
  const bool loaded_before = loader.is_loaded();
  const bool present = loader.available();
  REQUIRE(loader.is_loaded() == loaded_before);
  REQUIRE(loader.available() == present);
}

TEST_CASE("LammpsLoader: require_loaded is consistent with is_loaded",
          "[lammps][loader]") {
  auto &loader = eonc::LammpsLoader::instance();
  if (loader.is_loaded()) {
    REQUIRE(loader.open_no_mpi != nullptr);
    REQUIRE(loader.close != nullptr);
    REQUIRE(loader.command != nullptr);
    REQUIRE(loader.file != nullptr);
    REQUIRE(loader.scatter_atoms != nullptr);
    REQUIRE(loader.extract_variable != nullptr);
    REQUIRE_NOTHROW(loader.require_loaded());
  } else if (loader.available()) {
    // On disk but not yet loaded: the lazy load must succeed and fill in
    // every entry point.
    REQUIRE_NOTHROW(loader.require_loaded());
    REQUIRE(loader.is_loaded());
    REQUIRE(loader.open_no_mpi != nullptr);
    REQUIRE(loader.last_error().empty());
  } else {
    // available() is cwd / LD_LIBRARY_PATH only. dlopen of the soname
    // may still succeed via ld.so.cache. A hard miss must name a class.
    try {
      loader.require_loaded();
      REQUIRE(loader.is_loaded());
      REQUIRE(loader.open_no_mpi != nullptr);
      REQUIRE(loader.last_error().empty());
    } catch (const std::runtime_error &err) {
      const std::string what{err.what()};
      REQUIRE_THAT(what, Catch::Matchers::ContainsSubstring("liblammps"));
      REQUIRE_THAT(what,
                   Catch::Matchers::ContainsSubstring(loader.last_error()));
      REQUIRE_FALSE(loader.last_error().empty());
      REQUIRE_FALSE(loader.is_loaded());
      REQUIRE(loader.open_no_mpi == nullptr);
    }
  }
}

TEST_CASE("LAMMPS worker reap ignores EINTR", "[lammps][worker]") {
  REQUIRE(eonc::lammpsWorkerReaped(7, 7, 0));
  REQUIRE_FALSE(eonc::lammpsWorkerReaped(-1, 7, EINTR));
  REQUIRE(eonc::lammpsWorkerReaped(-1, 7, ECHILD));
}

TEST_CASE("LAMMPS open args honor logging", "[lammps][logging]") {
  const auto quiet = eonc::lammpsOpenArgs(false, false, "");
  REQUIRE(quiet[3] == "-log");
  REQUIRE(quiet[4] == "none");
  const auto logged = eonc::lammpsOpenArgs(true, true, "screen.tmp");
  REQUIRE(logged[1] == "-echo");
  bool saw_log_file = false;
  for (const auto &arg : logged) {
    if (arg == "log.lammps") {
      saw_log_file = true;
    }
  }
  REQUIRE_FALSE(saw_log_file);
  bool saw_screen = false;
  for (const auto &arg : logged) {
    if (arg == "screen.tmp") {
      saw_screen = true;
    }
  }
  REQUIRE(saw_screen);
  REQUIRE(logged.back() == "omp");
}

TEST_CASE("LAMMPS screen cursor rewinds when the file is replaced",
          "[lammps][logging]") {
  REQUIRE(eonc::lammpsScreenCursor(40, 10, false) == 0);
  REQUIRE(eonc::lammpsScreenCursor(40, 80, true) == 0);
  REQUIRE(eonc::lammpsScreenCursor(40, 80, false) == 40);
  REQUIRE(eonc::lammpsScreenCursor(-1, 80, false) == 0);
}

#ifndef _WIN32
TEST_CASE("LAMMPSPot Morse pair returns a finite energy", "[lammps][force]") {
  namespace fs = std::filesystem;
  auto &loader = eonc::LammpsLoader::instance();
  if (!loader.is_loaded()) {
    try {
      loader.require_loaded();
    } catch (const std::runtime_error &err) {
      const std::string what{err.what()};
      REQUIRE(what.find("liblammps") != std::string::npos);
      return;
    }
  }
  REQUIRE(loader.is_loaded());

  struct Cwd {
    fs::path previous;
    explicit Cwd(const fs::path &next) : previous(fs::current_path()) {
      fs::current_path(next);
    }
    ~Cwd() {
      std::error_code ec;
      fs::current_path(previous, ec);
    }
  };

  const auto writeInput = [](const fs::path &path, bool realUnits) {
    std::ofstream out(path);
    REQUIRE(out.is_open());
    if (realUnits) {
      out << "#!units real\n";
    }
    out << "pair_style morse 9.5\n";
    out << "pair_coeff * * 0.7102 1.6047 2.897\n";
    out << "pair_modify shift yes\n";
  };

  const double positions[6] = {0.0, 0.0, 0.0, 2.897, 0.0, 0.0};
  const int numbers[2] = {78, 78};
  const double box[9] = {20.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 20.0};
  const double mask[6] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  const fs::path metalDir =
      fs::temp_directory_path() /
      ("eon-lmp-metal-" + std::to_string(static_cast<long long>(::getpid())));
  fs::remove_all(metalDir);
  fs::create_directories(metalDir);
  writeInput(metalDir / "in.lammps", false);

  eonc::Parameters params;
  eonc::ParametersLoadAccess::potential_options(params).LAMMPSLogging = true;
  double metalEnergy = 0.0;
  double metalForces[6] = {};
  double warmEnergy = 0.0;
  double warmForces[6] = {};
  {
    const Cwd cwd(metalDir);
    LAMMPSPot pot(params);
    pot.setFixedMask(2, mask);
    pot.force(2, positions, numbers, metalForces, &metalEnergy, nullptr, box);
    pot.force(2, positions, numbers, warmForces, &warmEnergy, nullptr, box);
    if (pot.computesStress()) {
      REQUIRE(pot.cauchyStress().allFinite());
    }
  }
  REQUIRE(std::isfinite(metalEnergy));
  REQUIRE(metalEnergy < 1.0e5);
  REQUIRE(std::isfinite(metalForces[0]));
  REQUIRE(warmEnergy == Catch::Approx(metalEnergy).margin(1e-6));

  const fs::path realDir =
      fs::temp_directory_path() /
      ("eon-lmp-real-" + std::to_string(static_cast<long long>(::getpid())));
  fs::remove_all(realDir);
  fs::create_directories(realDir);
  writeInput(realDir / "in.lammps", true);
  double realEnergy = 0.0;
  double realForces[6] = {};
  {
    const Cwd cwd(realDir);
    LAMMPSPot pot(params);
    pot.force(2, positions, numbers, realForces, &realEnergy, nullptr, box);
  }
  REQUIRE(std::isfinite(realEnergy));
  REQUIRE(realEnergy * 23.0609 == Catch::Approx(metalEnergy).margin(1e-4));
  fs::remove_all(metalDir);
  fs::remove_all(realDir);
}
#endif

} // namespace tests
