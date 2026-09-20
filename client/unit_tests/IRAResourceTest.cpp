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

/// Catch2 mock for the IRA compare path. Does not dlopen libira.

#include "eon/libs/IRA/IRAResource.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/IRACompare.h"

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

namespace {
void stub_match(int /*nat1*/, const int * /*typ1*/, const double * /*coords1*/,
                const int * /*cand1*/, int nat2, const int * /*typ2*/,
                const double * /*coords2*/, const int * /*cand2*/,
                double /*distThreshold*/, double **rotMat, double **trans,
                int **perm, double *hausdorffDist, int *ierr) {
  if (rotMat && *rotMat) {
    double *R = *rotMat;
    for (int i = 0; i < 9; ++i) {
      R[i] = 0.0;
    }
    R[0] = R[4] = R[8] = 1.0;
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
    *hausdorffDist = 0.0;
  }
  if (ierr) {
    *ierr = 0;
  }
}

class MockIRAResource : public eonc::IIRAResource {
public:
  void require_loaded() override {}
  [[nodiscard]] bool is_loaded() const noexcept override { return true; }
  [[nodiscard]] libira_match_fn get_match_fn() const override {
    return &stub_match;
  }
  [[nodiscard]] libira_cshda_pbc_fn get_cshda_pbc_fn() const override {
    return nullptr;
  }
  [[nodiscard]] libira_compute_all_fn get_compute_all_fn() const override {
    return nullptr;
  }
  [[nodiscard]] libira_get_nmax_fn get_get_nmax_fn() const override {
    return nullptr;
  }
};
} // namespace

TEST_CASE("IRACompare matchArrays with injected mock does not load singleton",
          "[ira][resource][inject]") {
  const bool singleton_loaded = eonc::IRAResource::instance().is_loaded();
  MockIRAResource mock;
  const int typ[3] = {1, 1, 1};
  const double pos[9] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0};
  auto result =
      eonc::IRACompare::matchArrays(3, typ, pos, 3, typ, pos, 1.0, mock);
  REQUIRE(result.error == 0);
  REQUIRE_THAT(result.hausdorffDistance,
               Catch::Matchers::WithinAbs(0.0, 1e-12));
  REQUIRE(result.permutation.size() == 3);
  REQUIRE(result.permutation[0] == 0);
  REQUIRE(result.permutation[1] == 1);
  REQUIRE(result.permutation[2] == 2);
  REQUIRE(eonc::IRAResource::instance().is_loaded() == singleton_loaded);
}

} // namespace tests
