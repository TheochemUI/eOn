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
#include "catch2/catch_amalgamated.hpp"
#include "eon/Eigen.h"
#include "eon/SafeMath.h"
#include "eon/fpe_handler.h"

#include <cfenv>
#include <cmath>
#include <limits>

// The FPE handler must not re-trap the same instruction forever. Clearing
// sticky flags alone re-executes the faulting op with trapping still enabled;
// the handler must demote the exception mask so execution proceeds with Inf.
// This is the defect that flooded multi-GB stderr streams and froze clients
// mid-saddle-search once a single divide-by-zero occurred.
//
// Also covers SafeMath guards that keep application code off the trap path.

TEST_CASE("enableFPE then divide-by-zero continues once (no re-trap storm)",
          "[fpe]") {
#if defined(_WIN32)
  SKIP("Windows SEH FPE path covered separately");
#else
  eonc::enableFPE();

  volatile double num = 1.0;
  volatile double den = 0.0;
  volatile double result = 0.0;
  // Must not hang or abort: handler reports once and masks further ZE traps.
  result = num / den;

  REQUIRE(std::isinf(result));
  // Second div-by-zero must not re-enter a storm either (mask sticks).
  volatile double again = num / den;
  REQUIRE(std::isinf(again));

  eonc::disableFPE();
  // Restore a clean environment for later tests in the same process.
  feclearexcept(FE_ALL_EXCEPT);
#endif
}

TEST_CASE("safe_div returns fallback on zero denom without trapping",
          "[fpe][safemath]") {
  eonc::enableFPE();
  double v = eonc::safemath::safe_div(1.0, 0.0, 42.0);
  REQUIRE(v == Catch::Approx(42.0));
  v = eonc::safemath::safe_div(6.0, 2.0, -1.0);
  REQUIRE(v == Catch::Approx(3.0));
  eonc::disableFPE();
  feclearexcept(FE_ALL_EXCEPT);
}

TEST_CASE("safe_normalized returns zero for null vector under FPE traps",
          "[fpe][safemath]") {
  eonc::enableFPE();
  Eigen::Vector3d z = Eigen::Vector3d::Zero();
  Eigen::Vector3d n = eonc::safemath::safe_normalized(z);
  REQUIRE(n.norm() == Catch::Approx(0.0));
  Eigen::Vector3d u = eonc::safemath::safe_normalized(Eigen::Vector3d(3, 0, 0));
  REQUIRE(u.norm() == Catch::Approx(1.0));
  eonc::disableFPE();
  feclearexcept(FE_ALL_EXCEPT);
}

// Models the LAMMPS worker contract: client main arms traps, the forked pot
// worker must demote them before any external force eval (PairEAM::compute
// can raise FE_DIVBYZERO). Without demotion the child re-stormed under the
// old continue handler.
TEST_CASE("disableFPE demotes traps after enableFPE (worker path)",
          "[fpe][lammps]") {
#if defined(_WIN32)
  SKIP("Windows SEH path uses _controlfp_s; covered by enable/disable pair");
#elif defined(__unix__)
  eonc::enableFPE();
#if defined(FE_DIVBYZERO)
  REQUIRE((fegetexcept() & FE_DIVBYZERO) != 0);
#endif
  eonc::disableFPE();
#if defined(FE_DIVBYZERO)
  REQUIRE((fegetexcept() & FE_DIVBYZERO) == 0);
  REQUIRE((fegetexcept() & FE_INVALID) == 0);
  REQUIRE((fegetexcept() & FE_OVERFLOW) == 0);
#endif
  // Soft IEEE: no SIGFPE, result is Inf.
  volatile double r = 1.0 / 0.0;
  REQUIRE(std::isinf(r));
  feclearexcept(FE_ALL_EXCEPT);
#else
  SKIP("fegetexcept only on glibc/unix");
#endif
}

TEST_CASE("FPEHandler::eat_fpe demotes traps for external pot scopes",
          "[fpe][lammps]") {
#if defined(_WIN32) || !defined(__unix__)
  SKIP("eat_fpe / feholdexcept path exercised on unix");
#else
  eonc::enableFPE();
  {
    eonc::FPEHandler fpeh;
    fpeh.eat_fpe();
    // Non-stop environment: div-by-zero must not trap.
    volatile double r = 1.0 / 0.0;
    REQUIRE(std::isinf(r));
    fpeh.restore_fpe();
  }
  eonc::disableFPE();
  feclearexcept(FE_ALL_EXCEPT);
#endif
}
