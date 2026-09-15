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
//
// Krylov / Ritz space is the PHVA mobile set (see Lanczos.cpp): default
// phva_atoms=All keeps all free atoms; an explicit list bounds the space
// to 3 N_mobile without changing free/fixed.

#include "eon/Davidson.h"

#include "eon/EonLogger.h"
#include "eon/HelperFunctions.h"
#include "eon/MobileAtoms.h"
#include "eon/Potential.h"
#include "eon/SafeMath.h"
#include <algorithm>

#include <cmath>
#include <memory>
#include <vector>

namespace eonc {

Davidson::Davidson(std::shared_ptr<Matter> matter, const Parameters &params,
                   std::shared_ptr<Potential> pot)
    : LowestEigenmode(pot, params) {
  lowestEv.resize(matter->numberOfAtoms(), 3);
  lowestEv.setZero();
  lowestEw = 0.0;
}

void Davidson::compute(std::shared_ptr<Matter> matter, AtomMatrix direction) {
  const VectorXi mobile =
      resolveMobileAtoms(matter.get(), params.davidson_options().phva_atoms);
  compute(std::move(matter), std::move(direction), mobile);
}

void Davidson::compute(std::shared_ptr<Matter> matter, AtomMatrix direction,
                       const VectorXi &mobileIn) {
  totalForceCalls = 0;
  statsRotations = 0;
  lowestEv.resize(matter->numberOfAtoms(), 3);
  lowestEv.setZero();

  const VectorXi mobile = resolveMobileAtoms(matter.get(), mobileIn);
  const int size = 3 * static_cast<int>(mobile.size());
  if (size == 0) {
    lowestEw = 0.0;
    return;
  }

  const long maxIter = std::max(1L, params.davidson_options().max_iterations);
  const double tol = params.davidson_options().tolerance;
  const double dr = params.main_options().finiteDifference;
  const bool useDiagPrec = params.davidson_options().diagonal_preconditioner;

  MatrixXd V(size, maxIter);
  MatrixXd HV(size, maxIter);
  V.setZero();
  HV.setZero();

  VectorXd r = packMobileRows(direction, mobile);
  double beta = r.norm();
  if (beta < eonc::safemath::eps) {
    lowestEw = 0.0;
    return;
  }
  r /= beta;

  auto tmpMatter = std::make_unique<Matter>(*matter);
  const long forceCallsStart = tmpMatter->getForceCalls();
  const AtomMatrix pos0 = matter->getPositions();
  const VectorXd force0 = mobileForces(tmpMatter.get(), mobile);

  auto applyH = [&](const VectorXd &v) -> VectorXd {
    AtomMatrix pos = pos0;
    unpackMobileRows(packMobileRows(pos0, mobile) + dr * v, mobile, pos);
    tmpMatter->setPositions(pos);
    return -(mobileForces(tmpMatter.get(), mobile) - force0) / dr;
  };

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

    MatrixXd G = V.leftCols(subspace).transpose() * HV.leftCols(subspace);
    G = 0.5 * (G + G.transpose());

    Eigen::SelfAdjointEigenSolver<MatrixXd> es(G);
    ew = es.eigenvalues()(0);
    VectorXd y = es.eigenvectors().col(0);
    evEst = V.leftCols(subspace) * y;
    evEst.normalize();

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
                   "angle={:7.3f} dim={:3d} iter={:3d} n_mobile={}",
                   ew, ewAbsRelErr, residNorm, statsAngle, subspace, iter,
                   mobile.size());

    if (ewAbsRelErr < tol && residNorm < tol * (std::fabs(ew) + 1.0)) {
      QUILL_LOG_INFO(log, "[Davidson] Tolerance reached: {}", tol);
      break;
    }
    if (subspace >= maxIter) {
      QUILL_LOG_ERROR(log, "[Davidson] Max subspace dimension");
      break;
    }

    VectorXd t(size);
    if (useDiagPrec) {
      for (int k = 0; k < size; ++k) {
        const double denom = diagH(k) - ew;
        t(k) = resid(k) / (std::fabs(denom) > 1e-12 ? denom : 1e-12);
      }
    } else {
      t = resid;
    }

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
  lowestEv.setZero();
  if (evEst.size() == size) {
    unpackMobileRows(evEst, mobile, lowestEv);
  }
}

double Davidson::getEigenvalue() { return lowestEw; }

AtomMatrix Davidson::getEigenvector() { return lowestEv; }

} // namespace eonc
