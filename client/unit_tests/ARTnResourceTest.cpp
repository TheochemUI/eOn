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

/// Catch2 mock for the ARTn saddle path. Does not dlopen libartn.

#include "eon/libs/ARTn/ARTnResource.h"
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/ARTnSaddleSearch.h"
#include "eon/Matter.h"
#include "eon/Parameters.h"
#include "eon/Potential.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <memory>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

namespace {
int g_nat = 0;

int stub_create() { return 0; }
void stub_destroy() {}
int stub_set_param(const char *const /*name*/, const int /*crank*/,
                   const int * /*csize*/, const void * /*cval*/) {
  return 0;
}
void stub_setup(const int nat, bool *cerr) {
  g_nat = nat;
  if (cerr) {
    *cerr = false;
  }
}
void stub_step(const int nat, const double /*etot*/, double *const /*force*/,
               int const * /*ityp*/, double *const /*pos*/,
               const double * /*box*/, const int * /*if_pos*/,
               double * /*displ_vec*/, bool *lconv) {
  g_nat = nat;
  if (lconv) {
    *lconv = true;
  }
}
int stub_get_data(const char *name, void **cval) {
  if (name == nullptr || cval == nullptr) {
    return 1;
  }
  if (std::strcmp(name, "has_error") == 0) {
    auto *flag = static_cast<bool *>(std::malloc(sizeof(bool)));
    if (!flag) {
      return 1;
    }
    *flag = false;
    *cval = flag;
    return 0;
  }
  if (std::strcmp(name, "has_sad") == 0) {
    auto *flag = static_cast<bool *>(std::malloc(sizeof(bool)));
    if (!flag) {
      return 1;
    }
    *flag = true;
    *cval = flag;
    return 0;
  }
  if (std::strcmp(name, "tau_sad") == 0) {
    const std::size_t n = static_cast<std::size_t>(3 * g_nat);
    auto *tau = static_cast<double *>(std::malloc(n * sizeof(double)));
    if (!tau) {
      return 1;
    }
    for (std::size_t i = 0; i < n; ++i) {
      tau[i] = 0.0;
    }
    *cval = tau;
    return 0;
  }
  if (std::strcmp(name, "eigval_sad") == 0) {
    auto *eig = static_cast<double *>(std::malloc(sizeof(double)));
    if (!eig) {
      return 1;
    }
    *eig = -1.0;
    *cval = eig;
    return 0;
  }
  if (std::strcmp(name, "eigen_sad") == 0) {
    const std::size_t n = static_cast<std::size_t>(3 * g_nat);
    auto *evec = static_cast<double *>(std::malloc(n * sizeof(double)));
    if (!evec) {
      return 1;
    }
    for (std::size_t i = 0; i < n; ++i) {
      evec[i] = 0.0;
    }
    if (n > 0) {
      evec[0] = 1.0;
    }
    *cval = evec;
    return 0;
  }
  return 1;
}

int stub_get_data_no_tau(const char *name, void **cval) {
  if (name != nullptr && std::strcmp(name, "tau_sad") == 0) {
    if (cval) {
      *cval = nullptr;
    }
    return 1;
  }
  return stub_get_data(name, cval);
}

class MockARTnResource : public eonc::IARTnResource {
public:
  void require_loaded() override {}
  [[nodiscard]] bool is_loaded() const noexcept override { return true; }
  [[nodiscard]] artn_create_fn get_create_fn() const override {
    return &stub_create;
  }
  [[nodiscard]] setup_artn_fn get_setup_fn() const override {
    return &stub_setup;
  }
  [[nodiscard]] artn_fn get_artn_fn() const override { return nullptr; }
  [[nodiscard]] artn_destroy_fn get_destroy_fn() const override {
    return &stub_destroy;
  }
  [[nodiscard]] set_param_fn get_set_param_fn() const override {
    return &stub_set_param;
  }
  [[nodiscard]] get_param_fn get_get_param_fn() const override {
    return nullptr;
  }
  [[nodiscard]] get_runparam_fn get_get_runparam_fn() const override {
    return nullptr;
  }
  [[nodiscard]] get_data_fn get_get_data_fn() const override {
    return &stub_get_data;
  }
  [[nodiscard]] print_caller_fn get_print_caller_fn() const override {
    return nullptr;
  }
  [[nodiscard]] artn_step_fn get_artn_step_fn() const override {
    return &stub_step;
  }
  [[nodiscard]] get_error_fn get_get_error_fn() const override {
    return nullptr;
  }
};

class MockARTnResourceNoTau : public MockARTnResource {
public:
  [[nodiscard]] get_data_fn get_get_data_fn() const override {
    return &stub_get_data_no_tau;
  }
};
} // namespace

TEST_CASE("ARTnSaddleSearch run with injected mock does not load singleton",
          "[artn][resource][inject]") {
  const bool singleton_loaded = eonc::ARTnResource::instance().is_loaded();

  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->resize(2);
  VectorXi nrs(2);
  nrs << 1, 1;
  matter->setAtomicNrs(nrs);
  AtomMatrix pos = AtomMatrix::Zero(2, 3);
  pos(1, 0) = 1.5;
  matter->setPositions(pos);
  Matrix3d cell = Matrix3d::Identity() * 20.0;
  matter->setCell(cell);
  AtomMatrix mode = AtomMatrix::Zero(2, 3);
  mode(0, 0) = 1.0;

  MockARTnResource mock;
  auto search = std::make_unique<ARTnSaddleSearch>(matter, pot, mode, params);
  REQUIRE(search->run(mock) == ARTnSaddleSearch::STATUS_GOOD);
  REQUIRE(search->getForceCalls() > 0);
  REQUIRE_THAT(search->getEigenvalue(),
               Catch::Matchers::WithinAbs(-1.0, 1e-12));
  REQUIRE(eonc::ARTnResource::instance().is_loaded() == singleton_loaded);
}

TEST_CASE("ARTnSaddleSearch run refuses success without tau_sad",
          "[artn][resource][inject]") {
  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto matter = std::make_shared<Matter>(pot, params);
  matter->resize(2);
  VectorXi nrs(2);
  nrs << 1, 1;
  matter->setAtomicNrs(nrs);
  AtomMatrix pos = AtomMatrix::Zero(2, 3);
  pos(1, 0) = 1.5;
  matter->setPositions(pos);
  Matrix3d cell = Matrix3d::Identity() * 20.0;
  matter->setCell(cell);
  AtomMatrix mode = AtomMatrix::Zero(2, 3);
  mode(0, 0) = 1.0;

  MockARTnResourceNoTau mock;
  auto search = std::make_unique<ARTnSaddleSearch>(matter, pot, mode, params);
  REQUIRE(search->run(mock) == ARTnSaddleSearch::STATUS_BAD_ARTN_ERROR);
  REQUIRE(matter->getPositions()(1, 0) == 1.5);
  REQUIRE(std::isnan(search->getEigenvalue()));
}

} // namespace tests
