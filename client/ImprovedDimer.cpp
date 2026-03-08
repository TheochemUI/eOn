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
// An implementation of Johannes Kästner and Paul Sherwood's improved dimer.
// An attempt to keep to the variable names in their 2008 paper has been made.

#include "ImprovedDimer.h"
#include "HelperFunctions.h"
#include "LowestEigenmode.h"
#include "SafeMath.h"
#include "eonExceptions.hpp"

#include "EonLogger.h"
using namespace std;
using namespace helper_functions;

const char ImprovedDimer::OPT_SD[] = "sd";
const char ImprovedDimer::OPT_CG[] = "cg";
const char ImprovedDimer::OPT_LBFGS[] = "lbfgs";

ImprovedDimer::ImprovedDimer(std::shared_ptr<Matter> matter,
                             const Parameters &params,
                             std::shared_ptr<Potential> pot)
    : LowestEigenmode(pot, params) {
  x0 = std::make_shared<Matter>(pot, params);
  x1 = std::make_shared<Matter>(pot, params);
  *x0 = *matter;
  *x1 = *matter;
  tau.resize(3 * matter->numberOfAtoms());
  tau.setZero();
  totalForceCalls = 0;

  if (params.dimer_options.opt_method == OPT_CG) {
    init_cg = true;
  }
  /* Logger initialized via class member */
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
  // Handle initial direction tracking
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
  // RONEB setup
  // Track the best (most negative) curvature and corresponding mode
  double bestNegativeCurvature = std::numeric_limits<double>::max();
  VectorXd bestTau = tau;
  VectorXd bestX0Positions; // Positions at best curvature
  VectorXd bestG0;          // Gradient at best curvature
  VectorXd bestG1;          // Gradient at x1 for best curvature

  // Reference mode tracking
  VectorXd referenceMode;
  // This represents the "NEB Tangent" restricted to the free atoms.
  if (hasFixedReference) {
    // Use the persistent NEB tangent
    referenceMode = fixedReferenceMode;
  } else {
    // Standard Dimer behavior: preventing mode switching within one force call
    referenceMode = tau;
  }
  // --- end
  *x0 = *matter;
  *x1 = *matter;
  VectorXd x0_r = x0->getPositionsV();
  bestX0Positions = x0_r;

  x1->setPositionsV(x0_r + params.main_options.finiteDifference * tau);
  // Quick check: If we stepped into a wall, flip the tangent immediately
  // This saves the optimizer from fighting a huge gradient in the first step
  if (x1->getPotentialEnergy() - x0->getPotentialEnergy() >
      10.0 * params.main_options.finiteDifference) {
    tau = -tau;
    x1->setPositionsV(x0_r + params.main_options.finiteDifference * tau);
    QUILL_LOG_DEBUG(
        log, "[IDimer] Initial tangent flipped due to high energy wall.");
  }

  if (params.dimer_options.opt_method == OPT_LBFGS) {
    s.clear();
    y.clear();
    rho.clear();
    init_lbfgs = true;
  }

  // other vectors
  VectorXd x1_rp;
  VectorXd x1_r;
  VectorXd tau_prime;
  VectorXd tau_Old;
  VectorXd g1_prime;

  double delta = params.main_options.finiteDifference;
  double phi_tol =
      helper_functions::pi * (params.dimer_options.converged_angle / 180.0);
  double phi_prime = 0.0;
  double phi_min = 0.0;

  statsRotations = 0;

  auto x1p = std::make_shared<Matter>(pot, params);

  // Melander, Laasonen, Jonsson, JCTC, 11(3), 1055–1062, 2015
  // http://doi.org/10.1021/ct501155k
  if (params.dimer_options.remove_rotation) {
    rotationRemove(AtomMatrix::Map(x0_r.data(), x0->numberOfAtoms(), 3), x1);
    x1_r = x1->getPositionsV();
    tau = x1_r - x0_r;
    tau.normalize();
    x1_r = x0_r + tau * delta;
  }

  // Calculate the gradients on x0 and x1, g0 and g1, respectively.
  VectorXd g0 = -x0->getForcesV();
  VectorXd g1 = -x1->getForcesV();

  bestG0 = g0;
  bestG1 = g1;

  positions.clear();
  gradients.clear();

  positions.push_back(x0->getPositionsV());
  positions.push_back(x1->getPositionsV());
  gradients.push_back(g0);
  gradients.push_back(g1);

  do // while we have not reached phi_tol or maximum rotations.
  {

    // Calculate the rotational force, F_R.
    F_R = -2.0 * (g1 - g0) + 2.0 * ((g1 - g0).dot(tau)) * tau;

    statsTorque = F_R.norm() / (params.main_options.finiteDifference * 2.0);

    // Determine the step direction, theta
    if (params.dimer_options.opt_method == OPT_SD) // steepest descent
    {
      theta = eonc::safemath::safe_normalized(F_R);
    } else if (params.dimer_options.opt_method == OPT_CG) // conjugate gradients
    {
      if (init_cg) {
        init_cg = false;
        gamma = 0.0;
      } else {
        a = fabs(F_R.dot(F_R_Old));
        b = F_R_Old.squaredNorm();
        if (a < 0.5 * b) {
          gamma = F_R.dot(F_R - F_R_Old) / b;
        } else {
          gamma = 0.0;
        }
      }

      if (gamma == 0.0) {
        theta = F_R;
      } else {
        // theta = F_R + thetaOld * F_R_Old.norm() * gamma;
        // xph: use thetaOld before normalized, not F_R_Old
        theta = F_R + thetaOld * gamma;
      }

      theta = theta - theta.dot(tau) * tau;
      // xph: use thetaOld before normalized
      thetaOld = theta;
      double theta_norm = theta.norm();
      if (theta_norm < eonc::safemath::eps) {
        // theta parallel to tau; fall back to steepest descent direction
        theta = F_R - F_R.dot(tau) * tau;
      }
      theta.normalize();

      F_R_Old = F_R;
    } else if (params.dimer_options.opt_method == OPT_LBFGS) // quasi-newton
    {
      if (init_lbfgs == false) {
        // xph: s0 should the difference between tau and tau_Old, which are
        // normalized vectors. VectorXd s0 = x1->getPositionsV() - rPrev;
        VectorXd s0 = tau - tau_Old;
        s.push_back(s0);
        // xph: rescale the force; or rescale s0 = s0*delta
        VectorXd y0 = (F_R_Old - F_R) / delta;
        y.push_back(y0);
        rho.push_back(eonc::safemath::safe_recip(s0.dot(y0), 0.0));
      } else {
        init_lbfgs = false;
      }

      // xph: 60 is better than 10. The default H0 in ASE is 70.
      double H0 = 1. / 60.;

      int loopmax = s.size();
      std::vector<double> a(loopmax);

      VectorXd q = -F_R;

      for (int i = loopmax - 1; i >= 0; i--) {
        a[i] = rho[i] * s[i].dot(q);
        q -= a[i] * y[i];
      }

      VectorXd z = H0 * q;

      for (int i = 0; i < loopmax; i++) {
        double b = rho[i] * y[i].dot(z);
        z += s[i] * (a[i] - b);
      }

      double vd = -eonc::safemath::safe_normalized(z).dot(
          eonc::safemath::safe_normalized(F_R));
      if (vd > 1.0)
        vd = 1.0;
      if (vd < -1.0)
        vd = -1.0;
      double angle =
          eonc::safemath::safe_acos(vd) * (180.0 / helper_functions::pi);

      // xph: larger angle is allowed to avoid frequent restart (was 70.0)
      if (angle > 87.0) {
        s.clear();
        y.clear();
        rho.clear();
        z = -F_R;
      }

      theta = -eonc::safemath::safe_normalized(z);
      // xph:  rethogonalize theta to tau
      theta = theta - theta.dot(tau) * tau;
      eonc::safemath::safe_normalize_inplace(theta);

      thetaOld = theta;
      F_R_Old = F_R;
      tau_Old = tau;
    }

    // Calculate the curvature along tau, C_tau.
    C_tau = (g1 - g0).dot(tau) / delta;

    // Track the best negative curvature and its corresponding eigenvector
    if (C_tau < bestNegativeCurvature) {
      bestNegativeCurvature = C_tau;
      bestTau = tau;
      bestX0Positions = x0->getPositionsV();
      bestG0 = g0;
      bestG1 = g1;
      foundNegativeCurvature = (C_tau < 0.0);
    }

    // Calculate a rough estimate (phi_prime) of the optimum rotation angle.
    double d_C_tau_d_phi = 2.0 * (g1 - g0).dot(theta) / delta;
    phi_prime = -0.5 * eonc::safemath::safe_atan_ratio(d_C_tau_d_phi,
                                                       2.0 * abs(C_tau), 0.0);
    statsAngle = phi_prime * (180.0 / helper_functions::pi);

    double alignment = std::abs(tau.dot(referenceMode));
    if (abs(phi_prime) > phi_tol) {
      double b1 = 0.5 * d_C_tau_d_phi;

      // Calculate g1_prime.
      x0_r = x0->getPositionsV();
      // xph: renormalize the new tangent after rotating phi_prime
      // x1_rp = x0_r + (tau * cos(phi_prime) + theta * sin(phi_prime)) * delta;
      tau_prime = tau * cos(phi_prime) + theta * sin(phi_prime);
      tau_prime = eonc::safemath::safe_normalized(tau_prime);
      x1_rp = x0_r + tau_prime * delta;

      *x1p = *x1;
      x1p->setPositionsV(x1_rp);
      g1_prime = -x1p->getForcesV();

      // Update position and curvature histories
      positions.push_back(x1_rp);
      gradients.push_back(g1_prime);

      // Calculate C_tau_prime.
      // tau_prime = (x1_rp - x0_r) / (x1_rp - x0_r).norm(); //xph
      double C_tau_prime = (g1_prime - g0).dot(tau_prime) / delta;

      // Calculate phi_min.
      double a1 = eonc::safemath::safe_div(C_tau - C_tau_prime +
                                               b1 * sin(2.0 * phi_prime),
                                           1.0 - cos(2.0 * phi_prime), 0.0);
      double a0 = 2.0 * (C_tau - a1);
      phi_min = 0.5 * eonc::safemath::safe_atan_ratio(b1, a1, 0.0);

      // Determine the curvature for phi_min.
      double C_tau_min =
          0.5 * a0 + a1 * cos(2.0 * phi_min) + b1 * sin(2.0 * phi_min);

      // If the curvature is being maximized, push it over pi/2.
      if (C_tau_min > C_tau) {
        phi_min += helper_functions::pi * 0.5;
        C_tau_min =
            0.5 * a0 + a1 * cos(2.0 * phi_min) + b1 * sin(2.0 * phi_min);
      }

      // xph: for accurate LBFGS
      if (phi_min > helper_functions::pi * 0.5) {
        phi_min -= helper_functions::pi;
      }
      statsAngle = phi_min * (180.0 / helper_functions::pi);

      // Update x1, tau, and C_tau.
      // xph: normalize first
      tau = tau * cos(phi_min) + theta * sin(phi_min);
      tau = eonc::safemath::safe_normalized(tau);
      x1_r = x0_r + tau * delta;

      // Melander, Laasonen, Jonsson, JCTC, 11(3), 1055–1062, 2015
      // http://doi.org/10.1021/ct501155k
      if (params.dimer_options.remove_rotation) {
        x1->setPositionsV(x1_r);
        rotationRemove(AtomMatrix::Map(x0_r.data(), x0->numberOfAtoms(), 3),
                       x1);
        x1_r = x1->getPositionsV();
        tau = x1_r - x0_r;
        tau.normalize();
        x1_r = x0_r + tau * delta;
      }

      // x1_r = x0_r + (tau * cos(phi_min) + theta * sin(phi_min)) * delta;
      x1->setPositionsV(x1_r);
      // tau = (x1_r - x0_r) / (x1_r - x0_r).norm();

      C_tau = C_tau_min;

      // Calculate the new g1.
      double sin_phi_prime = sin(phi_prime);
      g1 = g1 * eonc::safemath::safe_div(sin(phi_prime - phi_min),
                                         sin_phi_prime, 0.0) +
           g1_prime *
               eonc::safemath::safe_div(sin(phi_min), sin_phi_prime, 0.0) +
           g0 * (1.0 - cos(phi_min) - sin(phi_min) * tan(phi_prime * 0.5));

      statsTorque = F_R.norm() / (2.0 * params.main_options.finiteDifference);
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
    if (alignment < params.neb_options.climbing_image.roneb.angle_tol &&
        params.neb_options.climbing_image.roneb.use_mmf) {
      QUILL_LOG_WARNING(
          log, "Terminating dimer due to lost mode (align {:.3f}).", alignment);

      rotationDidConverge = false;

      // Restore the best negative curvature state
      if (bestNegativeCurvature < 0.0) {
        C_tau = bestNegativeCurvature;
        tau = bestTau;

        // Restore x0 to the position where we had best curvature
        x0->setPositionsV(bestX0Positions);
        x1->setPositionsV(bestX0Positions + delta * bestTau);

        // Also update the input matter to reflect the restored position
        *matter = *x0;

        QUILL_LOG_DEBUG(log,
                        "Restored best negative curvature state:  C_tau={:.4f}",
                        C_tau);
        throw eonc::DimerModeRestoredException();
      } else {
        // Never found negative curvature - keep current (positive) value
        // to signal we're not at a saddle
        rotationDidConverge = false;
        QUILL_LOG_WARNING(
            log, "Never found negative curvature.  Final C_tau: {:.4f}", C_tau);
        throw eonc::DimerModeLostException();
      }
      break;
    }

  } while (abs(phi_prime) > abs(phi_tol) and abs(phi_min) > abs(phi_tol) and
           statsRotations < params.dimer_options.rotations_max);
}

double ImprovedDimer::getEigenvalue() { return C_tau; }

AtomMatrix ImprovedDimer::getEigenvector() {
  return AtomMatrix::Map(tau.data(), x0->numberOfAtoms(), 3);
}
