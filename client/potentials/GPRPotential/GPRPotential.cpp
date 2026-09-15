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

#include "eon/potentials/GPRPotential/GPRPotential.h"
#include "subprojects/gprdimer/gpr/auxiliary/AdditionalFunctionality.h"
#include "subprojects/gprdimer/structures/Structures.h"

#include <stdexcept>

GPRPotential::GPRPotential(const eonc::Parameters &p)
    : eonc::Potential(eonc::PotType::GPR, p) {
  gpr_model = nullptr;
}

void GPRPotential::registerGPRObject(
    gpr::GaussianProcessRegression *_gpr_model) {
  gpr_model = _gpr_model;
}

void GPRPotential::initialize(void) {}

void GPRPotential::cleanMemory(void) {}

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void GPRPotential::force(long N, const double *R, const int *atomicNrs,
                         double *F, double *U, double *variance,
                         const double *box) {
  variance = nullptr;
  gpr::Observation observation;

  // Copy R points. Note, R should correspond to the moving atoms only.
  observation.R.resize(1, N * 3);
  for (int i = 0; i < N; i++) {
    observation.R.set(i, {R[3 * i], R[3 * i + 1], R[3 * i + 2]});
  }

  // Note, the following functions should be called before calling for
  // gpr_model->calculatePotential() gpr_model->decomposeCovarianceMatrix(R,
  // ind) - takes covariance matrix and vector of repetitive indices
  // gpr_model->calculateMeanPrediction() - takes a vector of combined energy
  // and force gpr_model->calculatePosteriorMeanPrediction() - no arguments
  gpr_model->calculatePotential(observation);

  for (int i = 0; i < N; i++) {
    F[3 * i] = observation.G[3 * i];
    F[3 * i + 1] = observation.G[3 * i + 1];
    F[3 * i + 2] = observation.G[3 * i + 2];
  }

  if (observation.E.size() < 1) {
    throw std::runtime_error("GPRPotential: empty energy from GPR model");
  }
  *U = observation.E[0];
}
