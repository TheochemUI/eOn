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
#pragma once

/// Strongly-typed status codes used across optimizers and saddle-search
/// methods. Two enums:
///
///   StepResult    -- single-iteration outcome from
///                    Optimizer::step() / run() / line-search.
///   SaddleStatus  -- terminal status from any SaddleSearchMethod
///                    subclass (MinMode, ARTn, BasinHopping,
///                    DynamicsSaddleSearch, BGSD, ...).
///
/// Both keep stable underlying int values so downstream consumers
/// (results.dat parsers, the akmc.py driver, the integration tests
/// in JobIntegrationTest.cpp) that read the raw int continue to
/// work without re-interpretation. Use to_int() on the enum-class
/// values when you need to write the bare integer (e.g. format
/// strings, JSON, saved fixture files); use magic_enum::enum_name()
/// for the symbolic name in logs.

#include <string_view>

namespace eonc {

/// Single-iteration outcome from Optimizer::step() / Optimizer::run()
/// and any line-search loop layered on top. Pre-2.15 the same
/// information was carried as a bare int (-1 / 0 / 1) which made
/// every call site read like a flag check.
enum class StepResult : int {
  Failed = -1,      ///< Hard failure: NaN forces, internal error.
  NotConverged = 0, ///< Step taken but the convergence test is unmet.
  Converged = 1,    ///< Convergence achieved -- caller should stop.
};

[[nodiscard]] constexpr int to_int(StepResult r) noexcept {
  return static_cast<int>(r);
}

/// Saddle-search terminal status. Values mirror the historical
/// MinModeSaddleSearch::Status int constants byte-for-byte; the
/// only change is the type system. ARTn / BasinHopping /
/// DynamicsSaddleSearch / BGSD reuse the same set so callers get
/// one place to look up "what termination_reason did this job
/// emit". Numeric stability is enforced by explicit values.
///
/// New status codes append at the end and bump the
/// statusMessage() table.
enum class SaddleStatus : int {
  Good = 0,
  Init = 1,
  BadNoConvex = 2,
  BadHighEnergy = 3,
  BadMaxConcaveIterations = 4,
  BadMaxIterations = 5,
  BadNotConnected = 6,
  BadPrefactor = 7,
  BadHighBarrier = 8,
  BadMinima = 9,
  FailedPrefactor = 10,
  PotentialFailed = 11,
  NonnegativeAbort = 12,
  NonlocalAbort = 13,
  NegativeBarrier = 14,
  BadMdTrajectoryTooShort = 15,
  BadNoNegativeModeAtSaddle = 16,
  BadNoBarrier = 17,
  ZeromodeAbort = 18,
  OptimizerError = 19,
  DimerLostMode = 20,
  DimerRestoredBest = 21,
  /// pARTn's Fortran backend raised an error code outside the
  /// MinModeSaddleSearch::Status range. ARTn's pre-refactor
  /// constant was the magic 22.
  BadArtnError = 22,
};

[[nodiscard]] constexpr int to_int(SaddleStatus s) noexcept {
  return static_cast<int>(s);
}

/// Human-readable message for a SaddleStatus. Out-of-range integer
/// values returned by SaddleSearchMethod::getStatus() (e.g. from
/// future codes the binary doesn't know about) fall through to
/// "Unknown status".
[[nodiscard]] constexpr std::string_view statusMessage(SaddleStatus s) {
  switch (s) {
  case SaddleStatus::Good:
    return "Success";
  case SaddleStatus::Init:
    return "Initialized";
  case SaddleStatus::BadNoConvex:
    return "Initial displacement unable to reach convex region";
  case SaddleStatus::BadHighEnergy:
    return "Barrier too high";
  case SaddleStatus::BadMaxConcaveIterations:
    return "Too many iterations in concave region";
  case SaddleStatus::BadMaxIterations:
    return "Too many iterations";
  case SaddleStatus::BadNotConnected:
    return "Saddle is not connected to initial state";
  case SaddleStatus::BadPrefactor:
    return "Prefactors not within window";
  case SaddleStatus::BadHighBarrier:
    return "Energy barrier not within window";
  case SaddleStatus::BadMinima:
    return "Minimizations from saddle did not converge";
  case SaddleStatus::FailedPrefactor:
    return "Hessian calculation failed";
  case SaddleStatus::PotentialFailed:
    return "Potential evaluation failed";
  case SaddleStatus::NonnegativeAbort:
    return "Nonnegative initial mode, aborting";
  case SaddleStatus::NonlocalAbort:
    return "Nonlocal abort";
  case SaddleStatus::NegativeBarrier:
    return "Negative barrier detected";
  case SaddleStatus::BadMdTrajectoryTooShort:
    return "No reaction found during MD trajectory";
  case SaddleStatus::BadNoNegativeModeAtSaddle:
    return "Converged to stationary point with zero negative modes";
  case SaddleStatus::BadNoBarrier:
    return "No forward barrier found along minimized band";
  case SaddleStatus::ZeromodeAbort:
    return "Zero mode abort";
  case SaddleStatus::OptimizerError:
    return "Optimizer error";
  case SaddleStatus::DimerLostMode:
    return "Dimer lost mode";
  case SaddleStatus::DimerRestoredBest:
    return "Dimer restored best";
  case SaddleStatus::BadArtnError:
    return "ARTn library error";
  }
  return "Unknown status";
}

/// Backwards-compat overload for call sites that hold a raw int
/// (e.g. results.dat round-trip, integration tests that compare to
/// stored fixtures). Casts to SaddleStatus and dispatches.
[[nodiscard]] constexpr std::string_view statusMessage(int status) {
  return statusMessage(static_cast<SaddleStatus>(status));
}

} // namespace eonc
