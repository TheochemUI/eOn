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

/// Catch2 mock for the Metatomic C ABI host. Does not dlopen libmetatomic_pot.

#include "eon/potentials/Metatomic/MetatomicLoader.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/Parameters.h"
#include "eon/potentials/Metatomic/MetatomicDynPot.h"

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

namespace {
struct MockState {
  int require{0};
  int create{0};
  int destroy{0};
  int force{0};
  char handle{};
};

MockState g_mta;

EonMtaPot *stub_create(const EonMtaConfig *, char *, size_t) {
  ++g_mta.create;
  return reinterpret_cast<EonMtaPot *>(&g_mta.handle);
}

void stub_destroy(EonMtaPot *) { ++g_mta.destroy; }

int stub_force(EonMtaPot *, long nAtoms, const double *, const int *,
               double *forces, double *energy, double *variance,
               const double *) {
  ++g_mta.force;
  if (energy) {
    *energy = 0.0;
  }
  if (variance) {
    *variance = 0.0;
  }
  if (forces) {
    for (long i = 0; i < nAtoms * 3; ++i) {
      forces[i] = 0.0;
    }
  }
  return 0;
}

int stub_abi_version() { return EON_MTA_ABI_VERSION; }

class MockMetatomicLoader : public eonc::IMetatomicLoader {
public:
  MockMetatomicLoader() {
    create = &stub_create;
    destroy = &stub_destroy;
    force = &stub_force;
    abi_version = &stub_abi_version;
  }
  void require_loaded() override { ++g_mta.require; }
  [[nodiscard]] bool is_loaded() const noexcept override { return true; }
};
} // namespace

TEST_CASE("MetatomicDynPot constructs with injected loader mock",
          "[metatomic][loader][inject]") {
  const bool singleton_loaded = eonc::MetatomicLoader::instance().is_loaded();
  g_mta = MockState{};
  MockMetatomicLoader mock;
  eonc::Parameters params;
  double pos[3] = {0.0, 0.0, 0.0};
  int nrs[1] = {1};
  double forces[3] = {1.0, 1.0, 1.0};
  double energy = 1.0;
  double box[9] = {10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0};
  {
    eonc::MetatomicDynPot pot(params, mock);
    REQUIRE(g_mta.require == 1);
    REQUIRE(g_mta.create == 1);
    pot.force(1, pos, nrs, forces, &energy, nullptr, box);
    REQUIRE(g_mta.force == 1);
    REQUIRE_THAT(energy, Catch::Matchers::WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(forces[0], Catch::Matchers::WithinAbs(0.0, 1e-12));
  }
  REQUIRE(g_mta.destroy == 1);
  REQUIRE(eonc::MetatomicLoader::instance().is_loaded() == singleton_loaded);
}

TEST_CASE("MetatomicLoader: singleton returns consistent instance",
          "[metatomic][loader]") {
  auto &a = eonc::MetatomicLoader::instance();
  auto &b = eonc::MetatomicLoader::instance();
  REQUIRE(&a == &b);
}

} // namespace tests
