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

/// Catch2 mock for the PluginLoader path. Does not mutate the process
/// singleton search path or lib_present probe.

#include "eon/potentials/PluginLoader.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/potentials/RgpotAdapter/RgpotAdapter.h"
#include "rgpot/LennardJones/LJPot.hpp"

#include <stdexcept>
#include <vector>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

namespace {
class MockPluginLoader : public eonc::IPluginLoader {
public:
  int add_calls{0};
  void add_config_paths(const std::string &) override { ++add_calls; }
  [[nodiscard]] bool lib_present(const char *) const override { return false; }
  [[nodiscard]] const std::vector<std::string> &
  search_paths() const noexcept override {
    return paths_;
  }
  [[noreturn]] void throw_not_found(const char *, const char *) const override {
    throw std::runtime_error("mock plugin not found");
  }

protected:
  eonc::dynlib::Handle open_lib(const char *) override { return {}; }

private:
  std::vector<std::string> paths_{};
};
} // namespace

TEST_CASE("RgpotAdapter constructs with injected plugin loader mock",
          "[plugin][loader][inject]") {
  auto &singleton = eonc::PluginLoader::instance();
  const auto paths_before = singleton.search_paths();
  const bool present_before =
      singleton.lib_present("eon_no_such_potential_zzzz");

  MockPluginLoader mock;
  Parameters params;
  RgpotAdapter<rgpot::LJPot> pot(PotType::LJ, params, rgpot::LJConfig{}, mock);
  REQUIRE(mock.add_calls == 1);
  REQUIRE(singleton.search_paths() == paths_before);
  REQUIRE(singleton.lib_present("eon_no_such_potential_zzzz") ==
          present_before);
}

TEST_CASE("PluginLoader: singleton lib_present does not dlopen",
          "[plugin][loader]") {
  auto &loader = eonc::PluginLoader::instance();
  const auto paths_before = loader.search_paths();
  REQUIRE_FALSE(loader.lib_present("eon_no_such_potential_zzzz"));
  REQUIRE(loader.search_paths() == paths_before);
}

} // namespace tests
