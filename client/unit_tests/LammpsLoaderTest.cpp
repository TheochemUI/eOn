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

#include "potentials/LAMMPS/LammpsLoader.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

TEST_CASE("LammpsLoader: singleton returns consistent instance",
          "[lammps][loader]") {
  auto &a = eonc::LammpsLoader::instance();
  auto &b = eonc::LammpsLoader::instance();
  REQUIRE(&a == &b);
}

TEST_CASE("LammpsLoader: is_loaded reflects availability", "[lammps][loader]") {
  auto &loader = eonc::LammpsLoader::instance();
  // We cannot control whether LAMMPS is installed in the test environment,
  // but we can verify the interface is consistent.
  if (loader.is_loaded()) {
    // All required function pointers should be non-null
    REQUIRE(loader.open_no_mpi != nullptr);
    REQUIRE(loader.close != nullptr);
    REQUIRE(loader.command != nullptr);
    REQUIRE(loader.file != nullptr);
    REQUIRE(loader.scatter_atoms != nullptr);
    REQUIRE(loader.extract_variable != nullptr);
    // require_loaded should not throw
    REQUIRE_NOTHROW(loader.require_loaded());
  } else {
    // All function pointers should be null
    REQUIRE(loader.open_no_mpi == nullptr);
    // require_loaded should throw with a helpful message
    REQUIRE_THROWS_WITH(
        loader.require_loaded(),
        Catch::Matchers::ContainsSubstring("liblammps not found"));
  }
}

} // namespace tests
