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
#include "eon/NudgedElasticBand.h"
#include "eon/SolidStateNEB.h"

namespace eonc {
namespace {

long nebSegment(const NudgedElasticBand &neb) {
  return 3L * neb.atoms + (neb.solidState() ? 9L : 0L);
}

} // namespace

VectorXd NEBObjectiveFunction::getGradient(bool fdstep) {
  if (neb->movedAfterForceCall)
    neb->updateForces();
  const long seg = nebSegment(*neb);
  const long atomDof = 3L * neb->atoms;
  VectorXd gradV(seg * neb->numImages);
  for (long i = 1; i <= neb->numImages; i++) {
    // Negate in-place during copy to avoid a second pass over 40KB
    gradV.segment(seg * (i - 1), atomDof) =
        -VectorXd::Map(neb->projectedForce[i]->data(), atomDof);
    if (neb->solidState()) {
      const eonc::neb::CartesianStep step = eonc::neb::solidStateCartesianStep(
          *neb->path[i], *neb->projectedForce[i], neb->cellForce(i),
          neb->solidJacobian());
      gradV.segment(seg * (i - 1), atomDof) =
          -VectorXd::Map(step.positions.data(), atomDof);
      gradV.segment(seg * (i - 1) + atomDof, 9) =
          -VectorXd::Map(step.cell.data(), 9);
    }
  }
  return gradV;
}

double NEBObjectiveFunction::getEnergy() {
  // The band update evaluates dirty images as one batch; summing image
  // energies on a moved band would otherwise evaluate them one at a time.
  if (neb->movedAfterForceCall)
    neb->updateForces();
  double Energy{0};
  for (long i = 1; i <= neb->numImages; i++) {
    Energy += neb->path[i]->getPotentialEnergy();
  }
  return Energy;
}

void NEBObjectiveFunction::setPositions(const VectorXd &x) {
  neb->movedAfterForceCall = true;
  const long seg = nebSegment(*neb);
  const long atomDof = 3L * neb->atoms;
  for (long i = 1; i <= neb->numImages; i++) {
    const long offset = seg * (i - 1);
    if (neb->solidState()) {
      Matrix3d cell = Matrix3d::Map(x.segment(offset + atomDof, 9).data());
      cell(0, 1) = 0.0;
      cell(0, 2) = 0.0;
      cell(1, 2) = 0.0;
      neb->path[i]->setCell(cell);
    }
    neb->path[i]->setPositions(
        AtomMatrix::Map(x.segment(offset, atomDof).data(), neb->atoms, 3));
  }
}

VectorXd NEBObjectiveFunction::getPositions() {
  const long seg = nebSegment(*neb);
  const long atomDof = 3L * neb->atoms;
  VectorXd posV(seg * neb->numImages);
  for (long i = 1; i <= neb->numImages; i++) {
    const long offset = seg * (i - 1);
    posV.segment(offset, atomDof) =
        VectorXd::Map(neb->path[i]->getPositions().data(), atomDof);
    if (neb->solidState()) {
      posV.segment(offset + atomDof, 9) =
          VectorXd::Map(neb->path[i]->getCell().data(), 9);
    }
  }
  return posV;
}

int NEBObjectiveFunction::degreesOfFreedom() {
  return static_cast<int>(nebSegment(*neb) * neb->numImages);
}

bool NEBObjectiveFunction::isUncertain() {
  double maxMaxUnc = std::numeric_limits<double>::lowest();
  double currentMaxUnc{0};
  for (long idx = 0; idx <= neb->numImages + 1; idx++) {
    currentMaxUnc = neb->path[idx]->getEnergyVariance();
    if (currentMaxUnc > maxMaxUnc) {
      maxMaxUnc = currentMaxUnc;
    }
  }
  bool unc_conv{maxMaxUnc > params.gp_surrogate_options().uncertainty};
  if (unc_conv) {
    this->status = NudgedElasticBand::NEBStatus::MAX_UNCERTAINTY;
  }
  return unc_conv;
}

bool NEBObjectiveFunction::isConverged() {
  bool force_conv = getConvergence() < params.neb_options().force_tolerance;
  return force_conv;
}

double NEBObjectiveFunction::getConvergence() {
  return neb->convergenceForce();
}

VectorXd NEBObjectiveFunction::getMasses() const {
  long count = 0;
  for (long i = 1; i <= neb->numImages; ++i) {
    const Matter &image = *neb->path[static_cast<size_t>(i)];
    const AtomMatrix mask = image.getFree();
    for (long atom = 0; atom < image.numberOfAtoms(); ++atom) {
      if (mask.row(atom).sum() > 0.5) {
        ++count;
      }
    }
  }
  VectorXd masses(count);
  long written = 0;
  for (long i = 1; i <= neb->numImages; ++i) {
    const Matter &image = *neb->path[static_cast<size_t>(i)];
    const auto all = image.getMasses();
    const AtomMatrix mask = image.getFree();
    for (long atom = 0; atom < image.numberOfAtoms(); ++atom) {
      if (mask.row(atom).sum() > 0.5) {
        masses(written++) = all(atom);
      }
    }
  }
  return masses;
}

VectorXd NEBObjectiveFunction::difference(const VectorXd &a,
                                          const VectorXd &b) {
  const long seg = nebSegment(*neb);
  const long atomDof = 3L * neb->atoms;
  VectorXd pbcDiff(seg * neb->numImages);
  for (int i = 1; i <= neb->numImages; i++) {
    const int n = (i - 1) * static_cast<int>(seg);
    pbcDiff.segment(n, atomDof) =
        neb->path[i]->pbcV(a.segment(n, atomDof) - b.segment(n, atomDof));
    if (neb->solidState()) {
      pbcDiff.segment(n + atomDof, 9) =
          a.segment(n + atomDof, 9) - b.segment(n + atomDof, 9);
    }
  }
  return pbcDiff;
}

} // namespace eonc
