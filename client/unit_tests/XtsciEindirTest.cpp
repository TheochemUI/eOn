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

#include "eon/XtsciEindir.h"
#include "catch2/catch_amalgamated.hpp"

TEST_CASE("rgpot eindir descriptor matches the eOn objective contract",
          "[xtsci][eindir]") {
  const auto required = eonc::xtsci_eindir::eon_requirement();
  const auto rgpot = eonc::xtsci_eindir::rgpot_ev_angstrom();
  REQUIRE(eonc::xtsci_eindir::compatible(rgpot, required));
  auto hartree = rgpot;
  hartree.energy_unit = "hartree";
  REQUIRE_FALSE(eonc::xtsci_eindir::compatible(hartree, required));
  auto extra = rgpot;
  extra.operations |= 1ull << 3;
  REQUIRE(eonc::xtsci_eindir::compatible(extra, required));
  auto other_schema = rgpot;
  other_schema.schema_id = "other";
  REQUIRE_FALSE(eonc::xtsci_eindir::compatible(other_schema, required));
}

TEST_CASE("eindir ABI stamp must carry the gradient feature",
          "[xtsci][eindir]") {
  REQUIRE(eonc::xtsci_eindir::abi_accepts_gradient(
      1, 1, eonc::xtsci_eindir::kFeatureGradient));
  REQUIRE_FALSE(eonc::xtsci_eindir::abi_accepts_gradient(1, 1, 0));
  REQUIRE_FALSE(eonc::xtsci_eindir::abi_accepts_gradient(
      2, 1, eonc::xtsci_eindir::kFeatureGradient));
  auto wildcard = eonc::xtsci_eindir::eon_requirement();
  wildcard.producer_id.clear();
  wildcard.energy_sign = 0;
  REQUIRE(eonc::xtsci_eindir::compatible(
      eonc::xtsci_eindir::rgpot_ev_angstrom(), wildcard));
}
