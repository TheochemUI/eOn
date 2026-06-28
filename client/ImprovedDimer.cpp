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
// An implementation of Johannes Kaestner and Paul Sherwood's improved dimer.
// An attempt to keep to the variable names in their 2008 paper has been made.

#include "ImprovedDimer.h"
#include "DimerRotationDispatch.h"
#include "HelperFunctions.h"
#include "LowestEigenmode.h"
#include "SafeMath.h"
#include "eonExceptions.hpp"

#include <cmath>
#include <thread>

#include "EonLogger.h"
using namespace eonc::helpers;

const char ImprovedDimer::OPT_SD[] = "sd";
const char ImprovedDimer::OPT_CG[] = "cg";
const char ImprovedDimer::OPT_LBFGS[] = "lbfgs";

ImprovedDimer::ImprovedDimer(std::shared_ptr<Matter> matter,
                             const Parameters &params,
                             std::shared_ptr<Potential> pot)
    : LowestEigenmode(pot, params) {
  // Each dimer image gets its own potential for lock-free parallel evaluation
  auto x1Pot = (pot->needsPerImageInstance() && params.main_options.parallel)
                   ? eonc::helpers::makePotential(params)
                   : pot;
  x0 = std::make_shared<Matter>(pot, params);
  x1 = std::make_shared<Matter>(x1Pot, params);
  *x0 = *matter;
  *x1 = *matter;
  tau.resize(3 * matter->numberOfAtoms());
  tau.setZero();
  totalForceCalls = 0;

  if (params.dimer_options.opt_method == OPT_CG) {
    init_cg = true;
  }
}

void ImprovedDimer::setReferenceMode(const VectorXd &ref) {
  fixedReferenceMode = ref;
  if (fixedReferenceMode.norm() > 1e-10) {
    fixedReferenceMode.normalize();
  }
  hasFixedReference = true;
}

void ImprovedDimer::clearReferenceMode() { hasFixedReference = false; }

void ImprovedDimer::compute(std::shared_ptr<Matter> matter,
                            AtomMatrix initialDirectionAtomMatrix) {

  VectorXd initialDirection = VectorXd::Map(initialDirectionAtomMatrix.data(),
                                            3 * matter->numberOfAtoms());
  tau = initialDirection.array() * matter->getFreeV().array();
  rotationDidConverge = true;
  foundNegativeCurvature = false;
  if (tau.norm() > 1e-10) {
    tau.normalize();
  } else {
    // Fallback if the tangent was zero on free atoms (unlikely but safe)
    tau.setRandom();
    tau = tau.array() * matter->getFreeV().array();
    tau.normalize();
  }

  // Track the best (most negative) curvature and corresponding mode
  double bestNegativeCurvature = std::numeric_limits<double>::max();
  VectorXd bestTau = tau;
  VectorXd bestX0Positions;
  VectorXd bestG0, bestG1;

  // Reference mode tracking for OCINEB mode-switching prevention
  VectorXd referenceMode = hasFixedReference ? fixedReferenceMode : tau;

  *x0 = *matter;
  *x1 = *matter;
  VectorXd x0_r = x0->getPositionsV();
  bestX0Positions = x0_r;

  double delta = params.main_options.finiteDifference;
  x1->setPositionsV(x0_r + delta * tau);

  // If we stepped into a high-energy wall, flip the tangent immediately
  if (x1->getPotentialEnergy() - x0->getPotentialEnergy() > 10.0 * delta) {
    tau = -tau;
    x1->setPositionsV(x0_r + delta * tau);
    QUILL_LOG_DEBUG(
        log, "[IDimer] Initial tangent flipped due to high energy wall.");
  }

  // Optional: LOR / Lanczos / Davidson rotation backends (enum dispatch).
  if (auto alt = runAlternativeRotation(
          params.dimer_options.rotation_backend, matter, params, pot,
          AtomMatrix::Map(tau.data(), matter->numberOfAtoms(), 3),
          static_cast<quill::Logger *>(log))) {
    C_tau = alt->eigenvalue;
    tau = VectorXd::Map(alt->eigenvector.data(), 3 * matter->numberOfAtoms());
    totalForceCalls += alt->forceCalls;
    statsRotations = alt->rotations;
    tau = tau.array() * matter->getFreeV().array();
    if (tau.norm() > 1e-10) {
      tau.normalize();
    }
    x0_r = matter->getPositionsV();
    x0->setPositionsV(x0_r);
    x1->setPositionsV(x0_r + delta * tau);
    *matter = *x0;
    rotationDidConverge = true;
    foundNegativeCurvature = (C_tau < 0.0);
    return;
  }

  if (params.dimer_options.opt_method == OPT_LBFGS) {
    s.clear();
    y.clear();
    rho.clear();
    init_lbfgs = true;
  }

  VectorXd x1_rp, x1_r, tau_prime, tau_Old, g1_prime;
  double phi_tol =
      eonc::helpers::pi * (params.dimer_options.converged_angle / 180.0);
  double phi_prime = 0.0;
  double phi_min = 0.0;

  statsRotations = 0;

  // Use x1's potential for the trial rotation image (consistent with per-image)
  auto x1p = std::make_shared<Matter>(x1->getPotential(), params);

  // Melander, Laasonen, Jonsson, JCTC 11(3), 1055-1062, 2015
  if (params.dimer_options.remove_rotation) {
    rotationRemove(AtomMatrix::Map(x0_r.data(), x0->numberOfAtoms(), 3), x1);
    x1_r = x1->getPositionsV();
    tau = x1_r - x0_r;
    tau.normalize();
    x1_r = x0_r + tau * delta;
  }

  // Calculate gradients on x0 and x1.
  // Prefer batched evaluation when the potential supports it (single
  // model.forward() call for both replicas, e.g. MetatomicPotential on GPU).
  // Else fall back to thread-parallel when the potential is thread-safe or
  // wants per-image instances. Otherwise sequential.
  VectorXd g0, g1;
  bool canParallel =
      pot->isSharedInstanceThreadSafe() || pot->needsPerImageInstance();
  if (pot->supportsBatchEvaluation()) {
    long n = x0->numberOfAtoms();
    bool x0dirty = x0->needsForceUpdate();
    bool x1dirty = x1->needsForceUpdate();

    if (x0dirty && x1dirty) {
      auto nrs0 = x0->getAtomicNrs();
      auto nrs1 = x1->getAtomicNrs();
      auto box0 = x0->getCell();
      auto box1 = x1->getCell();
      const double *posVec[] = {x0->getPositions().data(),
                                x1->getPositions().data()};
      const int *nrsVec[] = {nrs0.data(), nrs1.data()};
      double *frcVec[] = {x0->forcesData(), x1->forcesData()};
      double energies[2], vars[2];
      const double *boxVec[] = {box0.data(), box1.data()};
      pot->forceBatch(2, n, posVec, nrsVec, frcVec, energies, vars, boxVec);
      x0->setComputedPotential(energies[0], vars[0]);
      x1->setComputedPotential(energies[1], vars[1]);
    } else if (x1dirty) {
      auto nrs = x1->getAtomicNrs();
      auto box = x1->getCell();
      const double *posVec[] = {x1->getPositions().data()};
      const int *nrsVec[] = {nrs.data()};
      double *frcVec[] = {x1->forcesData()};
      double energies[1], vars[1];
      const double *boxVec[] = {box.data()};
      pot->forceBatch(1, n, posVec, nrsVec, frcVec, energies, vars, boxVec);
      x1->setComputedPotential(energies[0], vars[0]);
    } else if (x0dirty) {
      x0->getForcesRaw(); // through computePotential
    }
    g0 = -x0->getForcesV();
    g1 = -x1->getForcesV();
  } else if (params.main_options.parallel && canParallel) {
    // std::thread instead of std::jthread (Apple Clang libc++). Guard so an
    // exception from the foreground call still joins t0 before rethrow.
    std::thread t0([&] { g0 = -x0->getForcesV(); });
    try {
      g1 = -x1->getForcesV();
    } catch (...) {
      if (t0.joinable())
        t0.join();
      throw;
    }
    t0.join();
  } else {
    g0 = -x0->getForcesV();
    g1 = -x1->getForcesV();
  }

  bestG0 = g0;
  bestG1 = g1;

  positions.clear();
  gradients.clear();
  positions.push_back(x0->getPositionsV());
  positions.push_back(x1->getPositionsV());
  gradients.push_back(g0);
  gradients.push_back(g1);

  do { // Rotation loop: converge phi or hit max rotations

    // Rotational force, F_R
    F_R = -2.0 * (g1 - g0) + 2.0 * ((g1 - g0).dot(tau)) * tau;
    statsTorque = F_R.norm() / (delta * 2.0);

    // Determine step direction theta via selected optimizer
    if (params.dimer_options.opt_method == OPT_SD) {
      theta = eonc::safemath::safe_normalized(F_R);

    } else if (params.dimer_options.opt_method == OPT_CG) {
      if (init_cg) {
        init_cg = false;
        gamma = 0.0;
      } else {
        a = std::abs(F_R.dot(F_R_Old));
        b = F_R_Old.squaredNorm();
        gamma = (a < 0.5 * b) ? F_R.dot(F_R - F_R_Old) / b : 0.0;
      }

      theta = (gamma == 0.0) ? F_R : F_R + thetaOld * gamma;
      theta -= theta.dot(tau) * tau;
      thetaOld = theta;
      if (theta.norm() < eonc::safemath::eps) {
        theta = F_R - F_R.dot(tau) * tau;
      }
      theta.normalize();
      F_R_Old = F_R;

    } else if (params.dimer_options.opt_method == OPT_LBFGS) {
      if (!init_lbfgs) {
        VectorXd s0 = tau - tau_Old;
        s.push_back(s0);
        VectorXd y0 = (F_R_Old - F_R) / delta;
        y.push_back(y0);
        rho.push_back(eonc::safemath::safe_recip(s0.dot(y0), 0.0));
      } else {
        init_lbfgs = false;
      }

      double H0 = 1.0 / 60.0;
      size_t loopmax = s.size();
      std::vector<double> alpha(loopmax);

      VectorXd q = -F_R;
      for (long i = static_cast<long>(loopmax) - 1; i >= 0; i--) {
        alpha[i] = rho[i] * s[i].dot(q);
        q -= alpha[i] * y[i];
      }
      VectorXd z = H0 * q;
      for (size_t i = 0; i < loopmax; i++) {
        double bv = rho[i] * y[i].dot(z);
        z += s[i] * (alpha[i] - bv);
      }

      double vd = std::clamp(-eonc::safemath::safe_normalized(z).dot(
                                 eonc::safemath::safe_normalized(F_R)),
                             -1.0, 1.0);
      double angle =
          eonc::safemath::safe_acos(vd) * (180.0 / eonc::helpers::pi);

      if (angle > 87.0) {
        s.clear();
        y.clear();
        rho.clear();
        z = -F_R;
      }

      theta = -eonc::safemath::safe_normalized(z);
      theta -= theta.dot(tau) * tau;
      eonc::safemath::safe_normalize_inplace(theta);

      thetaOld = theta;
      F_R_Old = F_R;
      tau_Old = tau;
    }

    // Curvature along tau
    C_tau = (g1 - g0).dot(tau) / delta;

    // Track best negative curvature for mode restoration
    if (C_tau < bestNegativeCurvature) {
      bestNegativeCurvature = C_tau;
      bestTau = tau;
      bestX0Positions = x0->getPositionsV();
      bestG0 = g0;
      bestG1 = g1;
      foundNegativeCurvature = (C_tau < 0.0);
    }

    // Estimate optimum rotation angle
    double d_C_tau_d_phi = 2.0 * (g1 - g0).dot(theta) / delta;
    phi_prime = -0.5 * eonc::safemath::safe_atan_ratio(
                           d_C_tau_d_phi, 2.0 * std::abs(C_tau), 0.0);
    statsAngle = phi_prime * (180.0 / eonc::helpers::pi);

    double alignment = std::abs(tau.dot(referenceMode));

    if (std::abs(phi_prime) > phi_tol) {
      double b1 = 0.5 * d_C_tau_d_phi;

      // Trial rotation to phi_prime
      x0_r = x0->getPositionsV();
      tau_prime = tau * std::cos(phi_prime) + theta * std::sin(phi_prime);
      tau_prime = eonc::safemath::safe_normalized(tau_prime);
      x1_rp = x0_r + tau_prime * delta;

      *x1p = *x1;
      x1p->setPositionsV(x1_rp);
      g1_prime = -x1p->getForcesV();

      positions.push_back(x1_rp);
      gradients.push_back(g1_prime);

      double C_tau_prime = (g1_prime - g0).dot(tau_prime) / delta;

      // Optimal rotation angle via Fourier interpolation
      double a1 = eonc::safemath::safe_div(
          C_tau - C_tau_prime + b1 * std::sin(2.0 * phi_prime),
          1.0 - std::cos(2.0 * phi_prime), 0.0);
      double a0 = 2.0 * (C_tau - a1);
      phi_min = 0.5 * eonc::safemath::safe_atan_ratio(b1, a1, 0.0);

      double C_tau_min = 0.5 * a0 + a1 * std::cos(2.0 * phi_min) +
                         b1 * std::sin(2.0 * phi_min);

      // If curvature is being maximized, push over pi/2
      if (C_tau_min > C_tau) {
        phi_min += eonc::helpers::pi * 0.5;
        C_tau_min = 0.5 * a0 + a1 * std::cos(2.0 * phi_min) +
                    b1 * std::sin(2.0 * phi_min);
      }

      // Keep phi_min in [-pi/2, pi/2] for accurate LBFGS
      if (phi_min > eonc::helpers::pi * 0.5) {
        phi_min -= eonc::helpers::pi;
      }
      statsAngle = phi_min * (180.0 / eonc::helpers::pi);

      // Apply optimal rotation
      tau = tau * std::cos(phi_min) + theta * std::sin(phi_min);
      tau = eonc::safemath::safe_normalized(tau);
      x1_r = x0_r + tau * delta;

      // Melander, Laasonen, Jonsson, JCTC 11(3), 1055-1062, 2015
      if (params.dimer_options.remove_rotation) {
        x1->setPositionsV(x1_r);
        rotationRemove(AtomMatrix::Map(x0_r.data(), x0->numberOfAtoms(), 3),
                       x1);
        x1_r = x1->getPositionsV();
        tau = x1_r - x0_r;
        tau.normalize();
        x1_r = x0_r + tau * delta;
      }

      x1->setPositionsV(x1_r);
      C_tau = C_tau_min;

      // Interpolate g1 at phi_min from g1 and g1_prime (saves one force call)
      double sin_pp = std::sin(phi_prime);
      g1 = g1 * eonc::safemath::safe_div(std::sin(phi_prime - phi_min), sin_pp,
                                         0.0) +
           g1_prime * eonc::safemath::safe_div(std::sin(phi_min), sin_pp, 0.0) +
           g0 * (1.0 - std::cos(phi_min) -
                 std::sin(phi_min) * std::tan(phi_prime * 0.5));

      statsTorque = F_R.norm() / (2.0 * delta);
      statsRotations += 1;
      QUILL_LOG_INFO(
          log,
          "[IDimerRot]  -----   ---------   ----------   ------------------   "
          "{:9.4f}   {:7.3f}   {:6.3f}   {:4}   {:5.3f}",
          C_tau, statsTorque, statsAngle, statsRotations, alignment);
    } else {
      QUILL_LOG_INFO(
          log,
          "[IDimerRot]  -----   ---------   ----------   ------------------   "
          "{:9.4f}   {:7.3f}   ------   ----   {:5.3f}",
          C_tau, F_R.norm() / delta, alignment);
    }

    // Check for mode loss (OCINEB dimer refinement)
    if (alignment < params.neb_options.climbing_image.ocineb.angle_tol &&
        params.neb_options.climbing_image.ocineb.use_mmf) {
      QUILL_LOG_WARNING(
          log, "Terminating dimer due to lost mode (align {:.3f}).", alignment);
      rotationDidConverge = false;

      if (bestNegativeCurvature < 0.0) {
        // Restore the best negative curvature state
        C_tau = bestNegativeCurvature;
        tau = bestTau;
        x0->setPositionsV(bestX0Positions);
        x1->setPositionsV(bestX0Positions + delta * bestTau);
        *matter = *x0;
        QUILL_LOG_DEBUG(
            log, "Restored best negative curvature state: C_tau={:.4f}", C_tau);
        throw eonc::DimerModeRestoredException();
      } else {
        QUILL_LOG_WARNING(
            log, "Never found negative curvature. Final C_tau: {:.4f}", C_tau);
        throw eonc::DimerModeLostException();
      }
    }

  } while (std::abs(phi_prime) > std::abs(phi_tol) &&
           std::abs(phi_min) > std::abs(phi_tol) &&
           statsRotations < params.dimer_options.rotations_max);
}

double ImprovedDimer::getEigenvalue() { return C_tau; }

AtomMatrix ImprovedDimer::getEigenvector() {
  return AtomMatrix::Map(tau.data(), x0->numberOfAtoms(), 3);
}
