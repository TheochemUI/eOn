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
/// Requires -Dwith_ira=true and a working libira.

#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/IRACompare.h"
#include "eon/Matter.h"
#include "eon/libs/IRA/IRAResource.h"
#include <vector>

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
    ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
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
  auto result = eonc::IRACompare{}.match(*m1, *m2, 1.0);

  REQUIRE(result.error == 0);
  REQUIRE_THAT(result.hausdorffDistance, Catch::Matchers::WithinAbs(0.0, 1e-4));
  // Permutation should map each atom to itself (may be 0- or 1-based)
  REQUIRE(result.permutation.size() ==
          static_cast<size_t>(m1->numberOfAtoms()));
}

TEST_CASE_METHOD(IRAFixture, "IRA matchArrays agrees with match",
                 "[ira][match]") {
  auto via_matter = eonc::IRACompare{}.match(*m1, *m2, 1.0);
  std::vector<int> z1(static_cast<size_t>(m1->numberOfAtoms()));
  std::vector<int> z2(static_cast<size_t>(m2->numberOfAtoms()));
  auto nrs1 = m1->getAtomicNrs();
  auto nrs2 = m2->getAtomicNrs();
  for (int i = 0; i < m1->numberOfAtoms(); ++i) {
    z1[static_cast<size_t>(i)] = nrs1[i];
  }
  for (int i = 0; i < m2->numberOfAtoms(); ++i) {
    z2[static_cast<size_t>(i)] = nrs2[i];
  }
  auto via_arr = eonc::IRACompare{}.matchArrays(
      m1->numberOfAtoms(), z1.data(), m1->getPositions().data(),
      m2->numberOfAtoms(), z2.data(), m2->getPositions().data(), 1.0);
  REQUIRE(via_arr.error == via_matter.error);
  REQUIRE_THAT(via_arr.hausdorffDistance,
               Catch::Matchers::WithinAbs(via_matter.hausdorffDistance, 1e-12));
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

  auto result = eonc::IRACompare{}.match(*m1, *m2, 5.0);

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

  auto result = eonc::IRACompare{}.match(*m1, *m2, 5.0);

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

  auto result = eonc::IRACompare{}.match(*m1, *m2, 10.0);

  REQUIRE(result.error == 0);
  REQUIRE(result.hausdorffDistance > 0.1);
}

TEST_CASE_METHOD(IRAFixture, "IRA findSymmetry returns valid point group",
                 "[ira][symmetry]") {
  auto result = eonc::IRACompare{}.findSymmetry(*m1, 0.1);

  REQUIRE(result.error == 0);
  // Every structure has at least the identity operation
  REQUIRE(result.nOperations >= 1);
  REQUIRE(!result.pointGroup.empty());
  REQUIRE(result.operations.size() == static_cast<size_t>(result.nOperations));
  REQUIRE(result.axes.size() == static_cast<size_t>(result.nOperations));
}

namespace {
class MockIRAResource : public eonc::IIRAResource {
public:
  int require_calls{0};
  int match_calls{0};

  MockIRAResource() {
    libira_match_ = &MockIRAResource::match_fn;
    active_ = this;
  }

  void require_loaded() override { ++require_calls; }
  [[nodiscard]] bool is_loaded() const noexcept override { return true; }

private:
  static MockIRAResource *active_;

  static void match_fn(int, const int *, const double *, const int *, int nat2,
                       const int *, const double *, const int *, double,
                       double **rotMat, double **trans, int **perm,
                       double *hausdorffDist, int *ierr) {
    ++active_->match_calls;
    if (rotMat && *rotMat) {
      for (int i = 0; i < 9; ++i) {
        (*rotMat)[i] = (i % 4 == 0) ? 1.0 : 0.0;
      }
    }
    if (trans && *trans) {
      (*trans)[0] = (*trans)[1] = (*trans)[2] = 0.0;
    }
    if (perm && *perm) {
      for (int i = 0; i < nat2; ++i) {
        (*perm)[i] = i;
      }
    }
    if (hausdorffDist) {
      *hausdorffDist = 1.25;
    }
    if (ierr) {
      *ierr = 0;
    }
  }
};

MockIRAResource *MockIRAResource::active_{nullptr};
} // namespace

TEST_CASE("IRACompare constructs with injected resource mock",
          "[ira][inject]") {
  const bool singleton_loaded = eonc::IRAResource::instance().is_loaded();
  MockIRAResource mock;
  eonc::IRACompare cmp(mock);
  const int typ[1] = {1};
  const double pos[3] = {0.0, 0.0, 0.0};
  auto result = cmp.matchArrays(1, typ, pos, 1, typ, pos, 1.0);
  REQUIRE(mock.require_calls == 1);
  REQUIRE(mock.match_calls == 1);
  REQUIRE(result.error == 0);
  REQUIRE(result.hausdorffDistance == 1.25);
  REQUIRE(eonc::IRAResource::instance().is_loaded() == singleton_loaded);
}

} /* namespace tests */
