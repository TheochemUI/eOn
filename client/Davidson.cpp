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
// Minimum-mode via Davidson subspace iteration with finite-difference
// Hessian-vector products (same H*v as Lanczos). Replaces dimer rotation
// constrained minimization when min_mode_method = davidson.

#include "Davidson.h"
#include "HelperFunctions.h"
#include "Potential.h"
#include "SafeMath.h"

#include "EonLogger.h"

#include <cmath>
#include <vector>

Davidson::Davidson(std::shared_ptr<Matter> matter, const Parameters &params,
                   std::shared_ptr<Potential> pot)
    : LowestEigenmode(pot, params) {
  lowestEv.resize(matter->numberOfAtoms(), 3);
  lowestEv.setZero();
  lowestEw = 0.0;
}

void Davidson::compute(std::shared_ptr<Matter> matter, AtomMatrix direction) {
  totalForceCalls = 0;
  statsRotations = 0;
  const int size = 3 * matter->numberOfFreeAtoms();
  const long maxIter = params.davidson_options.max_iterations;
  const double tol = params.davidson_options.tolerance;
  const double dr = params.main_options.finiteDifference;
  const bool useDiagPrec = params.davidson_options.diagonal_preconditioner;

  MatrixXd V(size, maxIter);
  MatrixXd HV(size, maxIter);
  V.setZero();
  HV.setZero();

  VectorXd r(size);
  int i, j;
  for (i = 0, j = 0; i < matter->numberOfAtoms(); i++) {
    if (!matter->getFixed(i)) {
      r.segment<3>(j) = direction.row(i);
      j += 3;
    }
  }

  double beta = r.norm();
  if (beta < eonc::safemath::eps) {
    lowestEw = 0.0;
    lowestEv.resize(matter->numberOfAtoms(), 3);
    lowestEv.setZero();
    return;
  }
  r /= beta;

  auto tmpMatter = std::make_unique<Matter>(*matter);
  const long forceCallsStart = tmpMatter->getForceCalls();
  const VectorXd force0 = tmpMatter->getForcesFreeV();

  // Apply H via one-sided FD of forces (same as Lanczos): H v ≈ -(F(x+dr
  // v)-F(x))/dr
  auto applyH = [&](const VectorXd &v) -> VectorXd {
    tmpMatter->setPositionsFreeV(matter->getPositionsFreeV() + dr * v);
    const VectorXd force1 = tmpMatter->getForcesFreeV();
    return -(force1 - force0) / dr;
  };

  // Optional cheap diagonal preconditioner: one FD probe along e_k is too
  // expensive, so approximate diag(H) from the first H*v0 components | (Hv)_i /
  // v_i | when |v_i| is appreciable; else 1.0.
  VectorXd diagH = VectorXd::Ones(size);

  V.col(0) = r;
  HV.col(0) = applyH(V.col(0));
  if (useDiagPrec) {
    for (int k = 0; k < size; ++k) {
      const double vk = V(k, 0);
      if (std::fabs(vk) > 1e-8) {
        diagH(k) = std::max(std::fabs(HV(k, 0) / vk), 1e-3);
      }
    }
  }

  double ew = 0.0, ewOld = 0.0;
  VectorXd evEst = V.col(0);
  VectorXd evOldEst = evEst;
  int subspace = 1;

  for (int iter = 0; iter < maxIter; ++iter) {
    statsRotations = iter;
    statsAngle = 0.0;

    // Rayleigh-Ritz on current subspace: G = V^T H V
    MatrixXd G = V.leftCols(subspace).transpose() * HV.leftCols(subspace);
    // Symmetrize numerical noise
    G = 0.5 * (G + G.transpose());

    Eigen::SelfAdjointEigenSolver<MatrixXd> es(G);
    ew = es.eigenvalues()(0);
    VectorXd y = es.eigenvectors().col(0);
    evEst = V.leftCols(subspace) * y;
    evEst.normalize();

    // Residual r = H x - lambda x  (using HV * y)
    VectorXd Hx = HV.leftCols(subspace) * y;
    VectorXd resid = Hx - ew * evEst;
    const double residNorm = resid.norm();
    const double ewAbsRelErr =
        (iter == 0) ? 1.0
                    : eonc::safemath::safe_div(std::fabs(ew - ewOld),
                                               std::fabs(ewOld), 1.0);
    ewOld = ew;
    statsTorque = std::max(ewAbsRelErr, residNorm / (std::fabs(ew) + 1e-12));
    statsAngle = eonc::safemath::safe_acos(std::fabs(evEst.dot(evOldEst))) *
                 (180 / eonc::helpers::pi);
    evOldEst = evEst;

    QUILL_LOG_INFO(log,
                   "[Davidson] ew={:10.6f} rel_err={:10.6f} |r|={:10.6f} "
                   "angle={:7.3f} dim={:3d} iter={:3d}",
                   ew, ewAbsRelErr, residNorm, statsAngle, subspace, iter);

    if (ewAbsRelErr < tol && residNorm < tol * (std::fabs(ew) + 1.0)) {
      QUILL_LOG_INFO(log, "[Davidson] Tolerance reached: {}", tol);
      break;
    }
    if (subspace >= maxIter) {
      QUILL_LOG_ERROR(log, "[Davidson] Max subspace dimension");
      break;
    }

    // Preconditioned residual correction (classical Davidson step)
    VectorXd t(size);
    if (useDiagPrec) {
      for (int k = 0; k < size; ++k) {
        const double denom = diagH(k) - ew;
        t(k) = resid(k) / (std::fabs(denom) > 1e-12 ? denom : 1e-12);
      }
    } else {
      t = resid;
    }

    // Orthogonalize against V and normalize
    for (int c = 0; c < subspace; ++c) {
      t -= V.col(c).dot(t) * V.col(c);
    }
    const double tnorm = t.norm();
    if (tnorm < 1e-14) {
      QUILL_LOG_ERROR(log, "[Davidson] Linear dependence in residual");
      break;
    }
    t /= tnorm;

    V.col(subspace) = t;
    HV.col(subspace) = applyH(t);
    if (useDiagPrec) {
      for (int k = 0; k < size; ++k) {
        const double vk = t(k);
        if (std::fabs(vk) > 1e-8) {
          const double d = std::fabs(HV(k, subspace) / vk);
          diagH(k) = std::max(diagH(k), std::max(d, 1e-3));
        }
      }
    }
    ++subspace;

    if (iter >= maxIter - 1) {
      QUILL_LOG_ERROR(log, "[Davidson] Max iterations");
      break;
    }
  }

  lowestEw = ew;
  totalForceCalls = tmpMatter->getForceCalls() - forceCallsStart;
  lowestEv.resize(matter->numberOfAtoms(), 3);
  lowestEv.setZero();
  for (i = 0, j = 0; i < matter->numberOfAtoms(); i++) {
    if (!matter->getFixed(i)) {
      lowestEv.row(i) = evEst.segment<3>(j);
      j += 3;
    }
  }
}

double Davidson::getEigenvalue() { return lowestEw; }

AtomMatrix Davidson::getEigenvector() { return lowestEv; }
