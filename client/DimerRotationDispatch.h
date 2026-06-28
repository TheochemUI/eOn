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

#include "BaseStructures.h"
#include "Davidson.h"
#include "EonLogger.h"
#include "LORRotation.h"
#include "Lanczos.h"
#include "Matter.h"
#include "Parameters.h"

#include <memory>
#include <optional>
#include <string_view>
#include <type_traits>

namespace eonc {

/// Result of a non-classical dimer rotation backend (mode estimate only).
struct DimerRotationResult {
  AtomMatrix eigenvector;
  double eigenvalue{0.0};
  long forceCalls{0};
  long rotations{0};
  /// LOR: true only when residual stop fired; Lanczos/Davidson: true (internal
  /// convergence criteria, not residual-exposed here).
  bool converged{true};
};

[[nodiscard]] inline constexpr bool
usesAlternativeRotation(DimerRotationBackend backend) noexcept {
  return backend != DimerRotationBackend::Classical;
}

namespace detail {

template <typename Solver>
[[nodiscard]] DimerRotationResult runRotationSolver(
    const std::shared_ptr<Matter> &matter, const Parameters &params,
    const std::shared_ptr<Potential> &pot, const AtomMatrix &initialDirection) {
  Solver solver(matter, params, pot);
  solver.compute(matter, initialDirection);
  DimerRotationResult out;
  out.eigenvector = solver.getEigenvector();
  out.eigenvalue = solver.getEigenvalue();
  out.forceCalls = solver.totalForceCalls;
  out.rotations = solver.statsRotations;
  if constexpr (std::is_same_v<Solver, LORRotation>) {
    out.converged = solver.convergedOnResidual;
  } else {
    out.converged = true;
  }
  return out;
}

} // namespace detail

/// Run Lanczos, Davidson, or LOR. Classical returns nullopt (stock rotation
/// loop).
[[nodiscard]] inline std::optional<DimerRotationResult> runAlternativeRotation(
    DimerRotationBackend backend, const std::shared_ptr<Matter> &matter,
    const Parameters &params, const std::shared_ptr<Potential> &pot,
    const AtomMatrix &initialDirection, quill::Logger *log = nullptr) {
  if (!usesAlternativeRotation(backend)) {
    return std::nullopt;
  }

  const std::string_view backendName = magic_enum::enum_name(backend);

  if (log) {
    QUILL_LOG_INFO(log,
                   "[DimerRot] rotation_backend={} (skipping classical "
                   "constrained rotation loop)",
                   backendName);
  }

  DimerRotationResult result;
  switch (backend) {
  case DimerRotationBackend::LOR:
    result = detail::runRotationSolver<LORRotation>(matter, params, pot,
                                                    initialDirection);
    break;
  case DimerRotationBackend::Lanczos:
    result = detail::runRotationSolver<Lanczos>(matter, params, pot,
                                                initialDirection);
    break;
  case DimerRotationBackend::Davidson:
    result = detail::runRotationSolver<Davidson>(matter, params, pot,
                                                 initialDirection);
    break;
  case DimerRotationBackend::Classical:
    return std::nullopt;
  }

  if (log) {
    QUILL_LOG_INFO(log,
                   "[DimerRot] alternative rotation done backend={} "
                   "eigenvalue={:.6f} force_calls={}",
                   backendName, result.eigenvalue, result.forceCalls);
  }
  return result;
}

} // namespace eonc
