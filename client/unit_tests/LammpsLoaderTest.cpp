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
#include "eon/potentials/PluginLoader.h"

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

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
  } else {
    // Absent from disk: the loader refuses by name and stays unloaded.
    REQUIRE_THROWS_WITH(
        loader.require_loaded(),
        Catch::Matchers::ContainsSubstring("liblammps not found"));
    REQUIRE_FALSE(loader.is_loaded());
    REQUIRE(loader.open_no_mpi == nullptr);
  }
}

TEST_CASE("PluginLoader: lib_present is a filesystem probe",
          "[plugin][loader]") {
  auto &loader = eonc::PluginLoader::instance();
  REQUIRE_FALSE(loader.lib_present("eon_no_such_potential_zzzz"));

  const auto tmp =
      std::filesystem::temp_directory_path() / "eon_plugin_probe_XXXXXX";
  std::error_code ec;
  std::filesystem::create_directories(tmp, ec);
  REQUIRE_FALSE(ec);
#ifdef _WIN32
  const std::string libname = "eon_probe_dummy.dll";
#else
#ifdef __APPLE__
  const std::string libname = "libeon_probe_dummy.dylib";
#else
  const std::string libname = "libeon_probe_dummy.so";
#endif
#endif
  {
    std::ofstream ofs(tmp / libname, std::ios::binary);
    ofs << "not a real library";
  }
  loader.add_config_paths(tmp.string());
  REQUIRE(loader.lib_present("eon_probe_dummy"));
  std::filesystem::remove_all(tmp, ec);
}

} // namespace tests
