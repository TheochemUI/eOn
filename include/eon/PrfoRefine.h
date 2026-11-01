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

#include "Matter.h"
#include "Parameters.h"

namespace eonc {

struct PrfoRefineResult {
  bool converged{false};
  /// The refined point had a larger max-atom force and was dropped.
  bool restored{false};
  long steps{0};
  long forceCalls{0};
  double maxForce{0.0};
  double energy{0.0};
};

/// Partitioned rational function optimization (P-RFO) of a first-order
/// saddle in the rgsaddle Sella order-1 session, configured by
/// [Saddle Search] prfo_refine_*. The session's Hessian model is the
/// identity with `mode` seeded at `curvature`. Only free atoms move, with
/// no rigid-body projection. `matter` ends at the last P-RFO point, or at
/// its starting point when that has the smaller max-atom force. A build
/// without rgsaddle, or against an application binary interface (ABI)
/// older than 1.15, throws std::runtime_error.
PrfoRefineResult prfoRefineSaddle(Matter &matter, const AtomMatrix &mode,
                                  double curvature, const Parameters &params);

/// Curvature along `mode` from one forward difference of the forces over
/// `dr` Angstrom. Costs one force call on a copy of `matter`.
double modeCurvature(const Matter &matter, const AtomMatrix &mode, double dr);

} // namespace eonc
