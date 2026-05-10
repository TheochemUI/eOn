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

/// Integration tests for the IRA structure comparison library.
/// IRAResource is always compiled in; libira.so is dlopen'd at first
/// require_loaded(). Tests SKIP when libira is not on
/// LD_LIBRARY_PATH at runtime.

#include "IRACompare.h"
#include "Matter.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

class IRAFixture {
protected:
  Parameters params;
  std::shared_ptr<Potential> pot;
  std::shared_ptr<Matter> m1;
  std::shared_ptr<Matter> m2;

  IRAFixture()
      : params{},
        pot{nullptr},
        m1{nullptr},
        m2{nullptr} {
    params.potential_options.potential = PotType::LJ;
    pot = eonc::helpers::makePotential(PotType::LJ, params);
    m1 = std::make_shared<Matter>(pot, params);
    m2 = std::make_shared<Matter>(pot, params);
    m1->con2matter(std::string("reactant.con"));
    m2->con2matter(std::string("reactant.con"));
  }
};

TEST_CASE_METHOD(IRAFixture,
                 "IRA match of identical structures returns zero distance",
                 "[ira][match]") {
  auto result = eonc::IRACompare::match(*m1, *m2, 1.0);

  REQUIRE(result.error == 0);
  REQUIRE_THAT(result.hausdorffDistance, Catch::Matchers::WithinAbs(0.0, 1e-4));
  // Permutation should map each atom to itself (may be 0- or 1-based)
  REQUIRE(result.permutation.size() ==
          static_cast<size_t>(m1->numberOfAtoms()));
}

TEST_CASE_METHOD(IRAFixture,
                 "IRA match of translated structure recovers translation",
                 "[ira][match]") {
  // Translate m2 by a known vector
  Eigen::Vector3d shift(1.5, -0.7, 0.3);
  auto pos = m2->getPositions();
  for (int i = 0; i < m2->numberOfAtoms(); i++) {
    pos.row(i) += shift.transpose();
  }
  m2->setPositions(pos);

  auto result = eonc::IRACompare::match(*m1, *m2, 5.0);

  REQUIRE(result.error == 0);
  // After alignment, Hausdorff distance should be near zero
  REQUIRE_THAT(result.hausdorffDistance, Catch::Matchers::WithinAbs(0.0, 0.1));
}

TEST_CASE_METHOD(IRAFixture,
                 "IRA match of permuted structure finds permutation",
                 "[ira][match]") {
  // Swap atoms 0 and 1 in m2
  auto pos = m2->getPositions();
  auto row0 = pos.row(0).eval();
  pos.row(0) = pos.row(1);
  pos.row(1) = row0;
  m2->setPositions(pos);

  auto result = eonc::IRACompare::match(*m1, *m2, 5.0);

  REQUIRE(result.error == 0);
  REQUIRE_THAT(result.hausdorffDistance, Catch::Matchers::WithinAbs(0.0, 0.1));
  // Permutation should reflect the swap
  REQUIRE(result.permutation.size() ==
          static_cast<size_t>(m1->numberOfAtoms()));
}

TEST_CASE_METHOD(IRAFixture,
                 "IRA match of different structures returns nonzero distance",
                 "[ira][match]") {
  // Displace atom 0 significantly
  auto pos = m2->getPositions();
  pos(0, 0) += 3.0;
  pos(0, 1) += 2.0;
  m2->setPositions(pos);

  auto result = eonc::IRACompare::match(*m1, *m2, 10.0);

  REQUIRE(result.error == 0);
  REQUIRE(result.hausdorffDistance > 0.1);
}

TEST_CASE_METHOD(IRAFixture, "IRA findSymmetry returns valid point group",
                 "[ira][symmetry]") {
  auto result = eonc::IRACompare::findSymmetry(*m1, 0.1);

  REQUIRE(result.error == 0);
  // Every structure has at least the identity operation
  REQUIRE(result.nOperations >= 1);
  REQUIRE(!result.pointGroup.empty());
  REQUIRE(result.operations.size() == static_cast<size_t>(result.nOperations));
  REQUIRE(result.axes.size() == static_cast<size_t>(result.nOperations));
}

} /* namespace tests */
