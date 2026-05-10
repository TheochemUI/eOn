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
#include "NudgedElasticBand.h"

namespace eonc {

VectorXd NEBObjectiveFunction::getGradient(bool fdstep) {
  if (neb->movedAfterForceCall)
    neb->updateForces();
  // long * int chain has to start in long so the implicit widening
  // doesn't trip bugprone-implicit-widening-of-multiplication-result
  // (atoms is `int` on the eOn side, seg / numImages are long).
  const long seg = static_cast<long>(3) * neb->atoms;
  VectorXd gradV(seg * neb->numImages);
  for (long i = 1; i <= neb->numImages; i++) {
    // Negate in-place during copy to avoid a second pass over 40KB
    gradV.segment(seg * (i - 1), seg) =
        -VectorXd::Map(neb->projectedForce[i]->data(), seg);
  }
  return gradV;
}

double NEBObjectiveFunction::getEnergy() {
  double Energy{0};
  for (long i = 1; i <= neb->numImages; i++) {
    Energy += neb->path[i]->getPotentialEnergy();
  }
  return Energy;
}

void NEBObjectiveFunction::setPositions(const VectorXd &x) {
  neb->movedAfterForceCall = true;
  const long seg = static_cast<long>(3) * neb->atoms;
  for (long i = 1; i <= neb->numImages; i++) {
    neb->path[i]->setPositions(
        AtomMatrix::Map(x.segment(seg * (i - 1), seg).data(), neb->atoms, 3));
  }
}

VectorXd NEBObjectiveFunction::getPositions() {
  VectorXd posV;
  const long seg = static_cast<long>(3) * neb->atoms;
  posV.resize(seg * neb->numImages);
  for (long i = 1; i <= neb->numImages; i++) {
    posV.segment(seg * (i - 1), seg) =
        VectorXd::Map(neb->path[i]->getPositions().data(), seg);
  }
  return posV;
}

int NEBObjectiveFunction::degreesOfFreedom() {
  return 3 * neb->numImages * neb->atoms;
}

bool NEBObjectiveFunction::isUncertain() {
  double maxMaxUnc = std::numeric_limits<double>::lowest();
  double currentMaxUnc{0};
  for (long idx = 0; idx <= neb->numImages; idx++) {
    currentMaxUnc = neb->path[idx]->getEnergyVariance();
    if (currentMaxUnc > maxMaxUnc) {
      maxMaxUnc = currentMaxUnc;
    }
  }
  bool unc_conv{maxMaxUnc > params.gp_surrogate_options.uncertainty};
  if (unc_conv) {
    this->status = NudgedElasticBand::NEBStatus::MAX_UNCERTAINTY;
  }
  return unc_conv;
}

bool NEBObjectiveFunction::isConverged() {
  bool force_conv = getConvergence() < params.neb_options.force_tolerance;
  return force_conv;
}

double NEBObjectiveFunction::getConvergence() {
  return neb->convergenceForce();
}

VectorXd NEBObjectiveFunction::difference(const VectorXd &a,
                                          const VectorXd &b) {
  VectorXd pbcDiff(3 * neb->numImages * neb->atoms);
  for (int i = 1; i <= neb->numImages; i++) {
    int n = (i - 1) * 3 * neb->atoms;
    int m = 3 * neb->atoms;
    pbcDiff.segment(n, m) =
        neb->path[i]->pbcV(a.segment(n, m) - b.segment(n, m));
  }
  return pbcDiff;
}

} // namespace eonc
