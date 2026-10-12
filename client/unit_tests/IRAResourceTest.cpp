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
#include "eon/Matter.h"
#include "eon/Parameters.h"
#include "eon/Potential.h"

#include <cstdlib>
#include <cstring>

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

void stub_cshda(int nat1, const int * /*typ1*/, const double * /*coords1*/,
                int /*nat2*/, const int * /*typ2*/, const double * /*coords2*/,
                const double * /*lat*/, double /*distThreshold*/, int **found,
                double **dists) {
  if (found != nullptr && *found != nullptr) {
    for (int i = 0; i < nat1; ++i) {
      (*found)[i] = i;
    }
  }
  if (dists != nullptr && *dists != nullptr) {
    for (int i = 0; i < nat1; ++i) {
      (*dists)[i] = 0.25;
    }
  }
}

int stub_nmax() { return 4; }

void stub_compute(int /*nat*/, const int * /*typ*/, const double * /*coords*/,
                  double /*threshold*/, int /*prescreenIh*/, int *n_mat,
                  double **mat, int ** /*perm*/, char **op, int ** /*n*/,
                  int ** /*p*/, double **ax, double **angle, double ** /*dH*/,
                  char **pg, int *n_prin_ax, double ** /*prin_ax*/, int *cerr) {
  if (n_mat != nullptr) {
    *n_mat = 1;
  }
  if (cerr != nullptr) {
    *cerr = 0;
  }
  if (n_prin_ax != nullptr) {
    *n_prin_ax = 1;
  }
  if (mat != nullptr) {
    auto *fresh = static_cast<double *>(std::malloc(9 * sizeof(double)));
    for (int i = 0; i < 9; ++i) {
      fresh[i] = 0.0;
    }
    fresh[0] = fresh[4] = fresh[8] = 1.0;
    *mat = fresh;
  }
  if (op != nullptr && *op != nullptr) {
    (*op)[0] = 'E';
    (*op)[1] = '\0';
  }
  if (ax != nullptr && *ax != nullptr) {
    (*ax)[0] = 0.0;
    (*ax)[1] = 0.0;
    (*ax)[2] = 1.0;
  }
  if (angle != nullptr && *angle != nullptr) {
    (*angle)[0] = 0.0;
  }
  if (pg != nullptr) {
    auto *fresh = static_cast<char *>(std::malloc(4));
    std::memcpy(fresh, "C1", 3);
    *pg = fresh;
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
    return &stub_cshda;
  }
  [[nodiscard]] libira_compute_all_fn get_compute_all_fn() const override {
    return &stub_compute;
  }
  [[nodiscard]] libira_get_nmax_fn get_get_nmax_fn() const override {
    return &stub_nmax;
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

TEST_CASE("IRACompare periodic assignment and symmetry use the injected resource",
          "[ira][resource][inject]") {
  eonc::Parameters params;
  auto pot = eonc::helpers::sharePotential(
      eonc::helpers::makePotential(eonc::PotType::LJ, params));
  eonc::Matter matter(pot, params);
  matter.resize(2);
  matter.setAtomicNr(0, 1);
  matter.setAtomicNr(1, 1);
  AtomMatrix positions(2, 3);
  positions.setZero();
  positions(1, 0) = 1.1;
  matter.setPositions(positions);

  MockIRAResource mock;
  const auto periodic =
      eonc::IRACompare::matchPBC(matter, matter, 0.5, mock);
  REQUIRE(periodic.error == 0);
  REQUIRE(periodic.permutation.size() == 2);
  REQUIRE(periodic.permutation[0] == 0);
  REQUIRE(periodic.permutation[1] == 1);
  REQUIRE(periodic.hausdorffDistance == Catch::Approx(0.25));

  const auto symmetry = eonc::IRACompare::findSymmetry(matter, 0.1, true, mock);
  REQUIRE(symmetry.error == 0);
  REQUIRE(symmetry.nOperations == 1);
  REQUIRE(symmetry.pointGroup == "C1");
  REQUIRE(symmetry.operations.size() == 1);
  REQUIRE(symmetry.axes.size() == 1);
  REQUIRE(symmetry.angles.size() == 1);

  eonc::Matter product(pot, params);
  product.resize(2);
  product.setAtomicNr(0, 1);
  product.setAtomicNr(1, 1);
  product.setPositions(positions);
  eonc::Matter reactant = matter;
  const auto aligned = eonc::IRACompare::alignReactantToProduct(
      reactant, product, 1.0, mock);
  REQUIRE(aligned.error == 0);
  REQUIRE(reactant.getPositions().isApprox(positions, 1e-12));
}

} // namespace tests
