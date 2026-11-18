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
#include "eon/Potential.h"
#include "eon/potentials/LAMMPS/LAMMPSPot.h"
#ifdef WITH_RGPOT
#include "eon/potentials/Rgpot/RgpotPot.h"
#endif

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__linux__)
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

TEST_CASE("LAMMPS fixed-mask updates ignore unchanged and all-free masks",
          "[lammps][mask][hotpath]") {
  const double freeMask[6] = {};
  const double fixedMask[6] = {0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
  std::vector<double> current;

  REQUIRE_FALSE(eonc::lammpsHasFixedCoordinates(2, freeMask));
  REQUIRE(eonc::lammpsHasFixedCoordinates(2, fixedMask));
  REQUIRE_FALSE(eonc::lammpsFixedMaskMatches(current, 2, freeMask));

  current.assign(fixedMask, fixedMask + 6);
  REQUIRE(eonc::lammpsFixedMaskMatches(current, 2, fixedMask));
  REQUIRE_FALSE(eonc::lammpsFixedMaskMatches(current, 2, freeMask));
}

#if defined(__linux__)
namespace {
std::vector<long> childProcesses() {
  std::ifstream in("/proc/self/task/" + std::to_string(getpid()) + "/children");
  std::vector<long> children;
  long pid = 0;
  while (in >> pid) {
    children.push_back(pid);
  }
  return children;
}
} // namespace

TEST_CASE("production LAMMPS stays in-process and matches the direct library",
          "[lammps][integration][inprocess]") {
  auto &loader = eonc::LammpsLoader::instance();
  try {
    loader.require_loaded();
  } catch (const std::runtime_error &) {
    SKIP("liblammps is not available");
  }

  const char *systems = std::getenv("EON_TEST_SYSTEMS_DIR");
  REQUIRE(systems != nullptr);
  const auto root =
      std::filesystem::weakly_canonical(std::filesystem::path(systems) /
                                        "../../../..");
  const auto table = root / "examples/CECAM_examples/advance_topic/FCC_GB/"
                            "potfiles/Cu01.eam.alloy";
  REQUIRE(std::filesystem::is_regular_file(table));

  const auto oldCwd = std::filesystem::current_path();
  const auto work = std::filesystem::temp_directory_path() /
                    ("eon-lammps-inprocess-" + std::to_string(getpid()));
  std::filesystem::create_directories(work);
  std::filesystem::current_path(work);
  {
    std::ofstream input("in.lammps");
    input << "pair_style eam/alloy\n"
          << "pair_coeff * * " << table.string() << " Cu\n";
  }

  eonc::Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LAMMPS;
  const auto before = childProcesses();
  LAMMPSPot production(params);
  REQUIRE(childProcesses() == before);
  LAMMPSPot direct(params, loader);
#ifdef WITH_RGPOT
  auto preferred = eonc::helpers::makePotential(params);
  auto *rgpot = dynamic_cast<RgpotPot *>(preferred.get());
  REQUIRE(rgpot != nullptr);
  REQUIRE(rgpot->backend() == "eam");

  eonc::Parameters rgpotLammpsParams;
  ParametersLoadAccess::potential_options(rgpotLammpsParams).potential =
      PotType::RGPOT;
  ParametersLoadAccess::rgpot_options(rgpotLammpsParams).backend = "lammps";
  auto rgpotLammps = eonc::helpers::makePotential(rgpotLammpsParams);
  REQUIRE(dynamic_cast<LAMMPSPot *>(rgpotLammps.get()) != nullptr);
  REQUIRE(childProcesses() == before);
#endif

  constexpr long n = 2;
  const double box[9] = {20.0, 0.0, 0.0, 0.0, 20.0,
                         0.0,  0.0, 0.0, 20.0};
  const double positions[6] = {5.0, 5.0, 5.0, 7.4, 5.0, 5.0};
  const int atomicNumbers[2] = {29, 29};
  double fp[6] = {}, fd[6] = {};
  double up = 0.0, ud = 0.0, variance = 0.0;
  production.force(n, positions, atomicNumbers, fp, &up, &variance, box);
  // An MPI-enabled liblammps may start one Open MPI singleton daemon on its
  // first library call. No subsequent force call may create another process.
  const auto afterLibraryInitialization = childProcesses();
  direct.force(n, positions, atomicNumbers, fd, &ud, &variance, box);

  REQUIRE(up == Catch::Approx(ud).margin(1.0e-8));
  for (int i = 0; i < 6; ++i) {
    REQUIRE(fp[i] == Catch::Approx(fd[i]).margin(1.0e-8));
  }
#ifdef WITH_RGPOT
  double fr[6] = {};
  double ur = 0.0;
  preferred->force(n, positions, atomicNumbers, fr, &ur, &variance, box);
  REQUIRE(ur == Catch::Approx(ud).margin(1.0e-8));
  for (int i = 0; i < 6; ++i) {
    REQUIRE(fr[i] == Catch::Approx(fd[i]).margin(1.0e-8));
  }
#endif
  REQUIRE(childProcesses() == afterLibraryInitialization);

  std::filesystem::current_path(oldCwd);
  std::filesystem::remove_all(work);
}
#endif

} // namespace tests
