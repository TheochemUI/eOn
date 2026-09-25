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

/// Tests for the metatomic runtime loader.
/// These tests verify an injected loader is used instead of the singleton.

#include "eon/potentials/Metatomic/MetatomicLoader.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Parameters.h"
#include "eon/potentials/Metatomic/MetatomicDynPot.h"

#include <string>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

namespace {
class MockMetatomicLoader : public eonc::IMetatomicLoader {
public:
  int require_calls{0};
  int force_calls{0};

  MockMetatomicLoader() {
    create = &MockMetatomicLoader::create_pot;
    destroy = &MockMetatomicLoader::destroy_pot;
    force = &MockMetatomicLoader::force_pot;
    active_ = this;
  }

  void require_loaded() override { ++require_calls; }
  [[nodiscard]] bool is_loaded() const noexcept override { return true; }
  bool try_load() override { return true; }

private:
  static MockMetatomicLoader *active_;

  static EonMtaPot *create_pot(const EonMtaConfig *, char *, size_t) {
    static int handle_token = 0;
    return reinterpret_cast<EonMtaPot *>(&handle_token);
  }
  static void destroy_pot(EonMtaPot *) {}
  static int force_pot(EonMtaPot *, long, const double *, const int *, double *,
                       double *, double *, const double *) {
    ++active_->force_calls;
    return 0;
  }
};

MockMetatomicLoader *MockMetatomicLoader::active_{nullptr};
} // namespace

TEST_CASE("MetatomicDynPot constructs with injected loader mock",
          "[metatomic][loader][inject]") {
  const bool singleton_loaded = eonc::MetatomicLoader::instance().is_loaded();
  MockMetatomicLoader mock;
  eonc::Parameters params;
  eonc::MetatomicDynPot pot(params, mock);
  REQUIRE(mock.require_calls == 1);
  double energy = 0.0;
  double pos[3] = {};
  int atomic = 1;
  double forces[3] = {};
  double box[9] = {};
  pot.force(1, pos, &atomic, forces, &energy, nullptr, box);
  REQUIRE(mock.force_calls == 1);
  REQUIRE(eonc::MetatomicLoader::instance().is_loaded() == singleton_loaded);
}

TEST_CASE("MetatomicLoader: singleton returns consistent instance",
          "[metatomic][loader]") {
  auto &a = eonc::MetatomicLoader::instance();
  auto &b = eonc::MetatomicLoader::instance();
  REQUIRE(&a == &b);
}

} // namespace tests
