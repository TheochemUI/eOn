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
#include "eon/potentials/LAMMPS/LAMMPSPot.h"

#include <stdexcept>
#include <string>

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

} // namespace tests
