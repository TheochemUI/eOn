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

#include "BondBoost.h"
#include "EonLogger.h"
#include "HelperFunctions.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

const char Hyperdynamics::NONE[] = "none";
const char Hyperdynamics::BOND_BOOST[] = "bond_boost";

BondBoost::BondBoost(Matter *matt, const Parameters &params)
    : matter{matt}, parameters{params} {
  nAtoms = matter->numberOfAtoms();
}

BondBoost::~BondBoost() = default;

void BondBoost::initialize() {
  nBBs = 0;
  nReg = 1;

  const std::string &balString = parameters.hyperdynamics_options.boost_atom_list;
  auto atoms = eonc::helpers::split_string_int(balString, ",");

  if (balString == "all" || atoms.empty()) {
    QUILL_LOG_DEBUG(log, "boost all atoms that are set free\n");
    nBAs = matter->numberOfFreeAtoms();
    nRAs = nAtoms - nBAs;
    BAList.resize(nBAs);
    RAList.resize(nRAs);
    long k = 0;
    for (long i = 0; i < nAtoms; i++) {
      if (!matter->getFixed(i)) {
        BAList[k++] = i;
      }
    }
  } else {
    QUILL_LOG_DEBUG(log, "boost the following selected atoms:");
    for (size_t i = 0; i < atoms.size(); i++) {
      QUILL_LOG_DEBUG(log, "{} ", atoms[i]);
    }
    QUILL_LOG_DEBUG(log, "\n");
    nBAs = static_cast<long>(atoms.size());
    nRAs = nAtoms - nBAs;
    BAList.resize(nBAs);
    RAList.resize(nRAs);
    for (long i = 0; i < nBAs; i++) {
      BAList[i] = atoms[i];
    }
  }

  // Build rest-atoms list (atoms not in BAList)
  long count = 0;
  for (long i = 0; i < nAtoms; i++) {
    bool isBoosted = std::any_of(BAList.begin(), BAList.end(),
                                 [i](long ba) { return ba == i; });
    if (!isBoosted) {
      RAList[count++] = i;
    }
  }
  if (count != nRAs) {
    QUILL_LOG_DEBUG(log,
                    "Error: nRestAtoms does not equal counted number!\n");
  }

  nTABs = nBAs * (nBAs - 1) / 2 + nBAs * nRAs;
  TABAList.resize(2 * nTABs);
  TABLList.setZero(nTABs, 1);
  QUILL_LOG_DEBUG(log, "BondBoost Used !\n");
}

double BondBoost::boost() {
  long RMDS = static_cast<long>(parameters.hyperdynamics_options.rmd_time /
                                parameters.dynamics_options.time_step);
  double biasPot = 0.0;

  if (nReg <= RMDS) {
    // Equilibration phase: accumulate average bond lengths
    Matrix<double, Eigen::Dynamic, 1> TABL_tmp = Rmdsteps();
    TABLList = TABLList + (1.0 / RMDS) * TABL_tmp;
    nReg++;
  } else {
    // Boost phase
    if (nReg == RMDS + 1) {
      nBBs = BondSelect();
    }
    Epsr_Q.resize(nBBs);
    CBBLList.setZero(nBBs, 1);
    biasPot = Booststeps();
    nReg++;
    Epsr_Q.clear();
  }
  return biasPot;
}

/// Compute bias potential and forces for bond-boost hyperdynamics.
/// Returns the bias potential energy contribution.
double BondBoost::Booststeps() {
  const double QRR = parameters.hyperdynamics_options.qrr;
  const double PRR = parameters.hyperdynamics_options.prr;
  const double DVMAX = parameters.hyperdynamics_options.dvmax;
  const double nBBsD = static_cast<double>(nBBs);

  AtomMatrix addForces(nBBs, 3);
  AtomMatrix TADF(nAtoms, 3);
  addForces.setZero();
  TADF.setZero();

  // Measure current bond lengths
  for (long i = 0; i < nBBs; i++) {
    CBBLList(i, 0) = matter->distance(BBAList[2 * i], BBAList[2 * i + 1]);
  }

  // Compute strain parameters and find maximum
  double epsrMax = 0.0;
  for (long i = 0; i < nBBs; i++) {
    Epsr_Q[i] = (CBBLList(i, 0) - EBBLList(i, 0)) / EBBLList(i, 0) / QRR;
    if (std::abs(Epsr_Q[i]) >= epsrMax) {
      epsrMax = std::abs(Epsr_Q[i]);
    }
  }

  // Envelope function A(eps_max)
  double A_eps = (1.0 - epsrMax * epsrMax) * (1.0 - epsrMax * epsrMax) /
                 (1.0 - PRR * PRR * epsrMax * epsrMax);

  // Sum of individual bias potentials
  double sumV = 0.0;
  if (epsrMax < 1.0) {
    for (long i = 0; i < nBBs; i++) {
      sumV += DVMAX * (1.0 - Epsr_Q[i] * Epsr_Q[i]) / nBBsD;
    }
  } else {
    A_eps = 0.0;
  }

  double boostFact = A_eps * sumV;

  // Compute bias forces per bond, accumulate on atoms
  for (long i = 0; i < nBBs; i++) {
    double dforce = 0.0;
    double fact1 = 2.0 * A_eps * DVMAX * Epsr_Q[i] / QRR /
                   EBBLList(i, 0) / nBBsD;

    if (std::abs(Epsr_Q[i]) < epsrMax) {
      dforce = fact1;
    } else {
      // Bond at maximum strain: additional envelope derivative
      double fTmp1 = 1.0 - PRR * PRR * Epsr_Q[i] * Epsr_Q[i];
      double fTmp2 = 1.0 - Epsr_Q[i] * Epsr_Q[i];
      double fact2 = 2.0 * fTmp2 * Epsr_Q[i] *
                     (2.0 * fTmp1 - PRR * PRR * fTmp2) / QRR /
                     EBBLList(i, 0) / fTmp1 / fTmp1;
      dforce = fact1 + sumV * fact2;
    }

    long a1 = BBAList[2 * i];
    long a2 = BBAList[2 * i + 1];
    double R = CBBLList(i, 0);

    for (int j = 0; j < 3; j++) {
      double rij = matter->pdistance(a1, a2, j);
      double fij = rij / R * dforce;
      addForces(i, j) = fij;
      TADF(a1, j) += fij;
      TADF(a2, j) -= fij;
    }
  }

  // Apply free-atom mask and set bias forces
  AtomMatrix biasForces = TADF.array() * matter->getFree().array();
  matter->setBiasForces(biasForces);
  return boostFact;
}

/// Measure all tagged-atom bond lengths for the equilibration phase.
Matrix<double, Eigen::Dynamic, 1> BondBoost::Rmdsteps() {
  Matrix<double, Eigen::Dynamic, 1> bondLengths(nTABs, 1);
  long count = 0;

  // Boost-atom pairs
  for (long i = 0; i < nBAs; i++) {
    for (long j = i + 1; j < nBAs; j++) {
      bondLengths(count, 0) = matter->distance(BAList[i], BAList[j]);
      TABAList[2 * count] = BAList[i];
      TABAList[2 * count + 1] = BAList[j];
      count++;
    }
  }

  // Boost-atom to rest-atom pairs
  for (long i = 0; i < nBAs; i++) {
    for (long j = 0; j < nRAs; j++) {
      bondLengths(count, 0) = matter->distance(BAList[i], RAList[j]);
      TABAList[2 * count] = BAList[i];
      TABAList[2 * count + 1] = RAList[j];
      count++;
    }
  }

  if (count != nTABs) {
    QUILL_LOG_DEBUG(
        log, "Total involved bond count does not match expected\n");
  }
  return bondLengths;
}

/// Select bonds within cutoff distance for boosting.
long BondBoost::BondSelect() {
  const double qCutoff = parameters.hyperdynamics_options.qcut;

  // Count bonds within cutoff
  long nSelected = 0;
  for (long i = 0; i < nTABs; i++) {
    if (TABLList(i, 0) <= qCutoff) {
      nSelected++;
    }
  }

  EBBLList.setZero(nSelected, 1);
  BBAList.resize(2 * nSelected);
  long count = 0;
  for (long i = 0; i < nTABs; i++) {
    if (TABLList(i, 0) <= qCutoff) {
      EBBLList(count, 0) = TABLList(i, 0);
      BBAList[2 * count] = TABAList[2 * i];
      BBAList[2 * count + 1] = TABAList[2 * i + 1];
      count++;
    }
  }
  return nSelected;
}
