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

#include <cstdlib>
#include <cstring>
#include <memory>
#include <mutex>
#include <thread>
#include <type_traits>
#include <vector>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

namespace {
int g_nat = 0;
int g_steps = 0;
int g_converge_after = 1;
std::mutex *g_mutex = nullptr;
int g_force_checks = 0;
int g_force_dropped = 0;
int g_destroy_calls = 0;
int g_nperp_calls = 0;
std::vector<int> g_nperp_vals;

void reset_artn_stubs() {
  g_nat = 0;
  g_steps = 0;
  g_converge_after = 1;
  g_mutex = nullptr;
  g_force_checks = 0;
  g_force_dropped = 0;
  g_destroy_calls = 0;
  g_nperp_calls = 0;
  g_nperp_vals.clear();
}

int stub_create() { return 0; }
void stub_destroy() { g_destroy_calls++; }
int stub_set_param(const char *const name, const int /*crank*/,
                   const int *csize, const void *cval) {
  if (name != nullptr && std::strcmp(name, "nperp_limitation") == 0 &&
      csize != nullptr && cval != nullptr && *csize > 0) {
    g_nperp_calls++;
    const auto *vals = static_cast<const int *>(cval);
    g_nperp_vals.assign(vals, vals + *csize);
  }
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
  g_steps++;
  if (lconv) {
    *lconv = g_steps >= g_converge_after;
  }
}
int stub_get_outptr(void **cval, int stored) {
  if (cval == nullptr) {
    return 1;
  }
  auto *value = static_cast<int *>(std::malloc(sizeof(int)));
  if (value == nullptr) {
    return 1;
  }
  *value = stored;
  *cval = value;
  return 0;
}

int stub_get_param(const char * /*name*/, void **cval) {
  return stub_get_outptr(cval, 7);
}

int stub_get_runparam(const char * /*name*/, void **cval) {
  return stub_get_outptr(cval, 11);
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
    return &stub_get_param;
  }
  [[nodiscard]] get_runparam_fn get_get_runparam_fn() const override {
    return &stub_get_runparam;
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
// Records whether library_mutex is free during force(). try_lock from the
// search thread is undefined if that thread already owns the mutex, so the
// probe runs on another thread.
class LockProbePotential : public eonc::Potential {
public:
  LockProbePotential()
      : eonc::Potential(PotType::LJ) {}

  void force(long nAtoms, const double * /*positions*/,
             const int * /*atomicNrs*/, double *forces, double *energy,
             double *variance, const double * /*box*/) override {
    for (long i = 0; i < nAtoms * 3; ++i) {
      forces[i] = 0.0;
    }
    if (energy) {
      *energy = 0.0;
    }
    if (variance) {
      *variance = 0.0;
    }
    if (g_mutex == nullptr) {
      return;
    }
    bool acquired = false;
    std::thread probe([&] {
      for (int attempt = 0; attempt < 32 && !acquired; ++attempt) {
        if (g_mutex->try_lock()) {
          acquired = true;
          g_mutex->unlock();
        }
      }
    });
    probe.join();
    g_force_checks += 1;
    if (acquired) {
      g_force_dropped += 1;
    }
  }
};

std::unique_ptr<ARTnSaddleSearch> make_artn_search(Parameters &params) {
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
  return std::make_unique<ARTnSaddleSearch>(matter, pot, mode, params);
}
} // namespace

TEST_CASE("ARTn get_param and get_runparam write through void**",
          "[artn][resource][get_param]") {
  using eonc::IARTnResource;
  static_assert(std::is_same_v<IARTnResource::get_param_fn,
                               int (*)(const char *, void **)>);
  static_assert(std::is_same_v<IARTnResource::get_runparam_fn,
                               int (*)(const char *, void **)>);

  MockARTnResource mock;
  void *param = nullptr;
  REQUIRE(mock.get_get_param_fn()("forc_thr", &param) == 0);
  REQUIRE(param != nullptr);
  REQUIRE(*static_cast<int *>(param) == 7);
  std::free(param);

  void *runparam = nullptr;
  REQUIRE(mock.get_get_runparam_fn()("PERP", &runparam) == 0);
  REQUIRE(runparam != nullptr);
  REQUIRE(*static_cast<int *>(runparam) == 11);
  std::free(runparam);
}

TEST_CASE("ARTnSaddleSearch run with injected mock does not load singleton",
          "[artn][resource][inject]") {
  reset_artn_stubs();
  const bool singleton_loaded = eonc::ARTnResource::instance().is_loaded();

  Parameters params;
  MockARTnResource mock;
  auto search = make_artn_search(params);
  REQUIRE(search->run(mock) == ARTnSaddleSearch::STATUS_GOOD);
  REQUIRE(search->getForceCalls() > 0);
  REQUIRE(g_nperp_calls == 0);
  REQUIRE(g_destroy_calls == 1);
  REQUIRE_THAT(search->getEigenvalue(),
               Catch::Matchers::WithinAbs(-1.0, 1e-12));
  REQUIRE(eonc::ARTnResource::instance().is_loaded() == singleton_loaded);
}

TEST_CASE("ARTnSaddleSearch holds library_mutex across force calls",
          "[artn][resource][lock]") {
  g_steps = 0;
  g_converge_after = 2;
  g_force_checks = 0;
  g_force_dropped = 0;
  g_mutex = nullptr;

  Parameters params;
  ParametersLoadAccess::potential_options(params).potential = PotType::LJ;
  auto pot = std::make_shared<LockProbePotential>();
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
  g_mutex = &mock.library_mutex;
  auto search = std::make_unique<ARTnSaddleSearch>(matter, pot, mode, params);
  const int status = search->run(mock);
  g_mutex = nullptr;

  REQUIRE(status == ARTnSaddleSearch::STATUS_GOOD);
  REQUIRE(g_steps >= 2);
  REQUIRE(g_force_checks >= 2);
  REQUIRE(g_force_dropped == 0);
}

TEST_CASE("ARTn nperp_limitation parse does not throw out of run",
          "[artn][resource][nperp]") {
  MockARTnResource mock;

  // "1," is not a blank token: getline stops at the trailing comma and
  // never yields an empty field. These are the tokens that throw.
  for (const char *bad : {",", " ", "1, ", "abc", "12abc", "2147483648"}) {
    reset_artn_stubs();
    Parameters params;
    ParametersLoadAccess::artn_options(params).nperp_limitation = bad;
    auto search = make_artn_search(params);
    int status = 0;
    REQUIRE_NOTHROW(status = search->run(mock));
    REQUIRE(status == ARTnSaddleSearch::STATUS_BAD_ARTN_ERROR);
    REQUIRE(search->getForceCalls() == 0);
    REQUIRE(g_destroy_calls == 1);
    REQUIRE(g_nperp_calls == 0);
  }

  reset_artn_stubs();
  Parameters params;
  ParametersLoadAccess::artn_options(params).nperp_limitation = "20, 30";
  auto search = make_artn_search(params);
  REQUIRE(search->run(mock) == ARTnSaddleSearch::STATUS_GOOD);
  REQUIRE(search->getForceCalls() > 0);
  REQUIRE(g_nperp_calls == 1);
  REQUIRE(g_nperp_vals == std::vector<int>{20, 30});
  REQUIRE(g_destroy_calls == 1);
}

} // namespace tests
