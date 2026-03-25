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
#include "NEBSpringForce.h"

namespace eonc::neb {

// --- UniformSpring ---

SpringResult
UniformSpring::compute(long i, const AtomMatrix &tangent, double distNext,
                       double distPrev, const AtomMatrix &posDiffNext,
                       const AtomMatrix &posDiffPrev,
                       const std::shared_ptr<Matter> &image) const {
  SpringResult result;
  result.forceSpringPar = ksp * (distNext - distPrev) * tangent;
  result.forceSpring = ksp * image->pbc((posDiffNext) - (posDiffPrev));
  return result;
}

// --- WeightedSpring ---

SpringResult WeightedSpring::compute(long i, const AtomMatrix &tangent,
                                     double distNext, double distPrev) const {
  double kspNext = springConstants[i];
  double kspPrev = springConstants[i - 1];

  SpringResult result;
  result.forceSpringPar =
      ((kspNext * distNext) - (kspPrev * distPrev)) * tangent;
  result.forceSpring = AtomMatrix::Zero(tangent.rows(), tangent.cols());
  return result;
}

// --- OnsagerMachlupSpring ---

SpringResult
OnsagerMachlupSpring::compute(long i, const AtomMatrix &tangent,
                              const AtomMatrix &posNext,
                              const AtomMatrix &posPrev, const AtomMatrix &pos,
                              const std::shared_ptr<Matter> &image) const {
  // Mandelli Eq. 13: k * ( R(i+1) + R(i-1) - 2R(i) + L(i+1) - L(i) )
  AtomMatrix diff =
      image->pbc(posNext + posPrev - 2.0 * pos + L_vecs[i + 1] - L_vecs[i]);
  AtomMatrix f_om_vec = base_k * diff;

  // Mandelli Eq. 15: Project onto tangent
  SpringResult result;
  result.forceSpringPar = matDot(f_om_vec, tangent) * tangent;
  result.forceSpring = AtomMatrix::Zero(tangent.rows(), tangent.cols());
  return result;
}

// --- Factory ---

SpringStrategy
buildSpringStrategy(const Parameters &params,
                    const std::vector<std::shared_ptr<Matter>> &path,
                    long numImages, int atoms, double maxEnergy, double E_ref) {

  if (params.neb_options.spring.om.enabled) {
    double base_k = params.neb_options.spring.constant;

    if (params.neb_options.spring.om.optimize_k) {
      double avgPotForce = 0.0;
      double avgPathCurvature = 0.0;
      int count = 0;
      for (long j = 1; j <= numImages; j++) {
        avgPotForce += path[j]->getForces().norm();
        AtomMatrix next = path[j + 1]->getPositions();
        AtomMatrix prev = path[j - 1]->getPositions();
        AtomMatrix curr = path[j]->getPositions();
        AtomMatrix curvVec = path[j]->pbc(next + prev - 2.0 * curr);
        avgPathCurvature += curvVec.norm();
        count++;
      }
      if (count > 0 && avgPathCurvature > 1e-6) {
        double scale = params.neb_options.spring.om.k_scale;
        base_k = scale * (avgPotForce / avgPathCurvature);
        base_k = std::max(base_k, params.neb_options.spring.om.k_min);
        base_k = std::min(base_k, params.neb_options.spring.om.k_max);
      }
    }

    // Pre-calculate L vectors for all images
    std::vector<AtomMatrix> L_vecs(numImages + 2);
    for (long j = 0; j <= numImages + 1; j++) {
      L_vecs[j].resize(atoms, 3);
      if (j == 0 || j == numImages + 1) {
        L_vecs[j].setZero();
      } else {
        const AtomMatrix &forces = path[j]->getForces();
        double alpha_k = eonc::safemath::safe_recip(2.0 * base_k, 0.0);
        for (int k = 0; k < atoms; k++) {
          L_vecs[j].row(k) = alpha_k * forces.row(k);
        }
      }
    }

    return OnsagerMachlupSpring{base_k, std::move(L_vecs)};

  } else if (params.neb_options.spring.weighting.enabled) {
    double k_l = params.neb_options.spring.weighting.k_min;
    double k_u = params.neb_options.spring.weighting.k_max;
    std::vector<double> springConstants(numImages + 2, k_l);

    double energyRange = maxEnergy - E_ref;
    if (energyRange < 1e-10) {
      std::fill(springConstants.begin(), springConstants.end(), k_l);
    } else {
      for (int idx = 1; idx <= numImages + 1; idx++) {
        double Ei = std::max(path[idx]->getPotentialEnergy(),
                             path[idx - 1]->getPotentialEnergy());
        if (Ei > E_ref) {
          double alpha_i = (maxEnergy - Ei) / energyRange;
          alpha_i = std::max(0.0, std::min(1.0, alpha_i));
          springConstants[idx - 1] = (1.0 - alpha_i) * k_u + alpha_i * k_l;
        } else {
          springConstants[idx - 1] = k_l;
        }
      }
    }

    return WeightedSpring{std::move(springConstants)};

  } else {
    return UniformSpring{params.neb_options.spring.constant};
  }
}

} // namespace eonc::neb
