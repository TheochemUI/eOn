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
#include "Eigen.h"
#include "EonLogger.h"

#include "LowestEigenmode.h"
#include "Matter.h"
#include "Parameters.h"

namespace eonc {

/// Davidson method for the lowest Hessian curvature mode (min-mode).
/// Uses the same finite-difference Hessian-vector products as Lanczos, but
/// expands a Ritz subspace with residual correction (optionally diagonal-
/// preconditioned) instead of dimer rotational constrained minimization.
/// See Olsen et al., J. Chem. Phys. 121, 9776 (2004) for FD context; Davidson
/// subspace expansion follows the classical residual-correction form.
class Davidson : public LowestEigenmode {

public:
  Davidson(std::shared_ptr<Matter> matter, const Parameters &params,
           std::shared_ptr<Potential> pot);
  ~Davidson() = default;
  void compute(std::shared_ptr<Matter> matter, AtomMatrix initialDirection);
  double getEigenvalue();
  AtomMatrix getEigenvector();

private:
  AtomMatrix lowestEv;
  double lowestEw;
  eonc::log::Scoped log;
};

} // namespace eonc

using eonc::Davidson;
