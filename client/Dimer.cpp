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
#include "Dimer.h"
#include "DimerRotationDispatch.h"
#include "HelperFunctions.h"
#include "SafeMath.h"

#include <cassert>
#include <cmath>
#include <thread>

using namespace eonc::helpers;

Dimer::Dimer(std::shared_ptr<Matter> matter, const Parameters &params,
             std::shared_ptr<Potential> pot)
    : LowestEigenmode(pot, params) {
  // Give matterDimer its own potential for parallel force evaluation
  auto dimerPot = (pot->needsPerImageInstance() && params.main_options.parallel)
                      ? eonc::helpers::makePotential(params)
                      : pot;
  matterCenter = std::make_shared<Matter>(pot, params);
  matterDimer = std::make_shared<Matter>(dimerPot, params);
  *matterCenter = *matter;
  *matterDimer = *matter;
  nAtoms = matter->numberOfAtoms();

  direction.resize(nAtoms, 3);
  rotationalPlane.resize(nAtoms, 3);
  direction.setZero();
  rotationalPlane.setZero();
  totalForceCalls = 0;
}

void Dimer::compute(std::shared_ptr<Matter> matter,
                    AtomMatrix initialDirection) {
  *matterCenter = *matter;

  eonc::safemath::safe_normalize_inplace(initialDirection);
  direction = initialDirection;

  // Optional: LOR / Lanczos / Davidson (enum dispatch; classical falls
  // through).
  if (auto alt = runAlternativeRotation(params.dimer_options.rotation_backend,
                                        matter, params, pot, direction,
                                        static_cast<quill::Logger *>(log))) {
    eigenvalue = alt->eigenvalue;
    direction = alt->eigenvector;
    totalForceCalls += alt->forceCalls;
    statsRotations = alt->rotations;
    eonc::safemath::safe_normalize_inplace(direction);
    *matterCenter = *matter;
    return;
  }

  long rotations = 0;
  long forceCallsCenter = matterCenter->getForceCalls();
  long forceCallsDimer = matterDimer->getForceCalls();
  double curvature = 0.0;
  double rotationAngle = 0.0;
  double torque = 0.0;

  AtomMatrix rotationalForce(nAtoms, 3);
  AtomMatrix rotationalForceOld(nAtoms, 3);
  AtomMatrix rotationalPlaneOld(nAtoms, 3);
  rotationalForce.setZero();
  rotationalForceOld.setZero();
  rotationalPlaneOld.setZero();

  statsAngle = 0;
  double lengthRotationalForceOld = 0.0;

  // Two force calls per rotation iteration
  bool doneRotating = false;
  while (!doneRotating) {
    curvature = calcRotationalForceReturnCurvature(rotationalForce);

    determineRotationalPlane(rotationalForce, rotationalForceOld,
                             rotationalPlaneOld, lengthRotationalForceOld);

    torque = rotationalForce.norm();
    assert(std::isnormal(torque));

    // Convergence: stop if torque is below threshold or max rotations reached
    if ((torque > params.dimer_options.torque_max &&
         rotations >= params.dimer_options.rotations_max) ||
        (torque < params.dimer_options.torque_max &&
         torque >= params.dimer_options.torque_min &&
         rotations >= params.dimer_options.rotations_min) ||
        (torque < params.dimer_options.torque_min)) {
      doneRotating = true;
    }

    double rotForce1 = matDot(rotationalForce, rotationalPlane);
    rotate(params.dimer_options.rotation_angle);

    if (!doneRotating) {
      curvature = calcRotationalForceReturnCurvature(rotationalForce);
      double rotForce2 = matDot(rotationalForce, rotationalPlane);

      double rotForceChange =
          (rotForce1 - rotForce2) / params.dimer_options.rotation_angle;
      double forceDimer = (rotForce1 + rotForce2) / 2.0;

      rotationAngle = eonc::safemath::safe_atan_ratio(2.0 * forceDimer,
                                                      rotForceChange, 0.0) /
                          2.0 -
                      params.dimer_options.rotation_angle / 2.0;

      if (rotForceChange < 0) {
        rotationAngle += eonc::helpers::pi / 2.0;
      }

      rotate(rotationAngle);
      rotationalPlaneOld = rotationalPlane;
      rotations++;
    }
    QUILL_LOG_DEBUG(log,
                    "[DimerRot]   -----   ---------   ----------------   "
                    "---------  {:9.3e}  {:9.3e}  {:9.3e}   ---------\n",
                    curvature, torque,
                    rotationAngle * (180.0 / eonc::helpers::pi));
  }

  statsTorque = torque;
  statsCurvature = curvature;
  direction.normalize();
  statsAngle = eonc::safemath::safe_acos(matDot(direction, initialDirection));
  statsAngle *= (180.0 / eonc::helpers::pi);
  statsRotations = rotations;
  eigenvalue = curvature;

  forceCallsCenter = matterCenter->getForceCalls() - forceCallsCenter;
  forceCallsDimer = matterDimer->getForceCalls() - forceCallsDimer;
  totalForceCalls += forceCallsCenter + forceCallsDimer;
}

double Dimer::getEigenvalue() { return eigenvalue; }

AtomMatrix Dimer::getEigenvector() { return direction; }

double Dimer::calcRotationalForceReturnCurvature(AtomMatrix &rotationalForce) {
  AtomMatrix posCenter = matterCenter->getPositions();

  // Displace to get dimer configuration A
  AtomMatrix posDimer =
      posCenter + direction * params.main_options.finiteDifference;

  // Optional rotation removal (Melander, Laasonen, Jonsson, JCTC 2015)
  if (params.dimer_options.remove_rotation) {
    matterDimer->setPositions(posDimer);
    rotationRemove(matterCenter, matterDimer);
    posDimer = matterDimer->getPositions();
    direction = posDimer - posCenter;
    direction.normalize();
    posDimer = posCenter + direction * params.main_options.finiteDifference;
  }

  // Obtain forces for dimer and center
  matterDimer->setPositions(posDimer);
  AtomMatrix forceA, forceCenter;
  if (pot->supportsBatchEvaluation()) {
    // Only batch systems that actually need recomputation.
    bool centerDirty = matterCenter->needsForceUpdate();
    bool dimerDirty = matterDimer->needsForceUpdate();

    if (centerDirty && dimerDirty) {
      // Both need eval -- batch together
      auto nrs0 = matterCenter->getAtomicNrs();
      auto nrs1 = matterDimer->getAtomicNrs();
      auto box0 = matterCenter->getCell();
      auto box1 = matterDimer->getCell();
      const double *posVec[] = {matterCenter->getPositions().data(),
                                posDimer.data()};
      const int *nrsVec[] = {nrs0.data(), nrs1.data()};
      double *frcVec[] = {matterCenter->forcesData(),
                          matterDimer->forcesData()};
      double energies[2], vars[2];
      const double *boxVec[] = {box0.data(), box1.data()};
      pot->forceBatch(2, nAtoms, posVec, nrsVec, frcVec, energies, vars,
                      boxVec);
      matterCenter->setComputedPotential(energies[0], vars[0]);
      matterDimer->setComputedPotential(energies[1], vars[1]);
    } else if (dimerDirty) {
      // Only dimer moved -- eval just dimer, center is cached
      auto nrs = matterDimer->getAtomicNrs();
      auto box = matterDimer->getCell();
      const double *posVec[] = {posDimer.data()};
      const int *nrsVec[] = {nrs.data()};
      double *frcVec[] = {matterDimer->forcesData()};
      double energies[1], vars[1];
      const double *boxVec[] = {box.data()};
      pot->forceBatch(1, nAtoms, posVec, nrsVec, frcVec, energies, vars,
                      boxVec);
      matterDimer->setComputedPotential(energies[0], vars[0]);
    } else if (centerDirty) {
      // Only center moved (rare)
      matterCenter->getForces(); // through computePotential
    }
    // else: both cached, nothing to do
    forceCenter = matterCenter->getForces();
    forceA = matterDimer->getForces();
  } else {
    forceA = matterDimer->getForces();
    forceCenter = matterCenter->getForces();
  }
  AtomMatrix forceB = 2.0 * forceCenter - forceA;

  double projA = matDot(direction, forceA);
  double projB = matDot(direction, forceB);

  // Remove force component parallel to dimer
  forceA = makeOrthogonal(forceA, direction);
  forceB = makeOrthogonal(forceB, direction);

  // Rotational force = orthogonal force difference
  rotationalForce =
      (forceA - forceB) / (2.0 * params.main_options.finiteDifference);

  // Curvature along the dimer
  return (projB - projA) / (2.0 * params.main_options.finiteDifference);
}

void Dimer::determineRotationalPlane(const AtomMatrix &rotationalForce,
                                     AtomMatrix &rotationalForceOld,
                                     const AtomMatrix &rotationalPlaneOld,
                                     double &lengthRotationalForceOld) {
  double gamma = 0.0;
  double a = std::abs(matDot(rotationalForce, rotationalForceOld));
  double b = rotationalForceOld.squaredNorm();
  if (a < 0.5 * b) {
    // Polak-Ribiere conjugate gradient direction
    gamma = (rotationalForce.array() *
             (rotationalForce - rotationalForceOld).array())
                .sum() /
            b;
  }

  // New rotational plane from current force and previous plane
  rotationalPlane =
      rotationalForce + rotationalPlaneOld * lengthRotationalForceOld * gamma;

  // Orthogonalize to dimer direction and normalize
  lengthRotationalForceOld = rotationalPlane.norm();
  rotationalPlane = makeOrthogonal(rotationalPlane, direction);
  eonc::safemath::safe_normalize_inplace(rotationalPlane);

  rotationalForceOld = rotationalForce;
}

void Dimer::rotate(double rotationAngle) {
  statsAngle += rotationAngle;

  double cosA = std::cos(rotationAngle);
  double sinA = std::sin(rotationAngle);

  direction = direction * cosA + rotationalPlane * sinA;
  rotationalPlane = rotationalPlane * cosA - direction * sinA;

  eonc::safemath::safe_normalize_inplace(direction);
  eonc::safemath::safe_normalize_inplace(rotationalPlane);

  // Remove component from rotationalPlane parallel to direction
  rotationalPlane = makeOrthogonal(rotationalPlane, direction);
  eonc::safemath::safe_normalize_inplace(rotationalPlane);
}
