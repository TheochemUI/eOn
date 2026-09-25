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
// This Lanczos algorithm is implemented as described in this paper:
// R. A. Olsen, G. J. Kroes, G. Henkelman, A. Arnaldsson, and H. Jónsson,
// Comparison of methods for finding saddle points without knowledge of the
// final states, J. Chem. Phys. 121, 9776-9792 (2004).
//
// Krylov space is the PHVA mobile set (Li & Jensen, Theor. Chem. Acc. 107,
// 211 (2002)): free/fixed is the optimizer mask; phva_atoms (or an explicit
// mobile list) selects which free atoms enter the 3 N_mobile Krylov space.
// Default phva_atoms=All keeps historical behavior (all free atoms).

#include "eon/Lanczos.h"
#include "eon/EonLogger.h"
#include "eon/FiniteDifference.h"
#include "eon/HelperFunctions.h"
#include "eon/MobileAtoms.h"
#include "eon/Potential.h"
#include "eon/SafeMath.h"

#include <algorithm>
#include <cmath>
#include <memory>

namespace eonc {

Lanczos::Lanczos(std::shared_ptr<Matter> matter, const Parameters &params,
                 std::shared_ptr<Potential> pot)
    : LowestEigenmode(pot, params),
      lowestEw{0.0} {
  lowestEv.resize(matter->numberOfAtoms(), 3);
  lowestEv.setZero();
}

void Lanczos::compute(std::shared_ptr<Matter> matter, AtomMatrix direction) {
  const VectorXi mobile =
      resolveMobileAtoms(matter.get(), params.lanczos_options().phva_atoms);
  compute(std::move(matter), std::move(direction), mobile);
}

void Lanczos::compute(std::shared_ptr<Matter> matter, AtomMatrix direction,
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

  const long maxIters = std::max(1L, params.lanczos_options().max_iterations);
  MatrixXd T(size, maxIters), Q(size, maxIters);
  T.setZero();
  VectorXd u(size), r = packMobileRows(direction, mobile);

  double alpha, beta = r.norm();
  if (beta < eonc::safemath::eps) {
    lowestEw = 0.0;
    return;
  }
  double ew = 0, ewOld = 0, ewAbsRelErr;
  const double dr = params.main_options().finiteDifference;
  const FdScheme requested = parseFdScheme(params.hessian_options().fd_scheme);
  // Fourth-order is the only scheme that changes the matrix-free product.
  // one_sided and central keep the forward difference used by min-mode.
  const FdScheme hvpScheme =
      requested == FdScheme::Fourth ? FdScheme::Fourth : FdScheme::OneSided;
  VectorXd evEst, evT, evOldEst;

  auto tmpMatter = std::make_unique<Matter>(*matter);
  const long forceCallsStart = tmpMatter->getForceCalls();
  const AtomMatrix pos0 = matter->getPositions();
  const VectorXd force0 = mobileForces(tmpMatter.get(), mobile);

  auto applyH = [&](const VectorXd &v) -> VectorXd {
    auto at = [&](double scale) -> VectorXd {
      AtomMatrix pos = pos0;
      unpackMobileRows(packMobileRows(pos0, mobile) + (scale * dr) * v, mobile,
                       pos);
      tmpMatter->setPositions(pos);
      return mobileForces(tmpMatter.get(), mobile);
    };
    return fdHessianVector(hvpScheme, dr, force0, at);
  };

  for (int i = 0; i < size; i++) {
    statsRotations = i;
    Q.col(i) = r / beta;

    u = applyH(Q.col(i));

    if (i == 0) {
      r = u;
    } else {
      r = u - beta * Q.col(i - 1);
    }
    alpha = Q.col(i).dot(r);
    r = r - alpha * Q.col(i);

    T(i, i) = alpha;
    if (i > 0) {
      T(i - 1, i) = beta;
      T(i, i - 1) = beta;
    }

    beta = r.norm();

    if (beta <= 1e-10 * std::fabs(alpha)) {
      if (i == 0) {
        ew = alpha;
        evEst = Q.col(0);
      }
      QUILL_LOG_ERROR(log, "[ILanczos] ERROR: linear dependence");
      break;
    }
    if (i >= 1) {
      Eigen::SelfAdjointEigenSolver<MatrixXd> es(T.block(0, 0, i + 1, i + 1));
      ew = es.eigenvalues()(0);
      evT = es.eigenvectors().col(0);
      ewAbsRelErr = eonc::safemath::safe_div(std::fabs(ew - ewOld),
                                             std::fabs(ewOld), 1.0);
      ewOld = ew;

      evEst = Q.block(0, 0, size, i + 1) * evT;
      evEst.normalize();
      statsAngle = eonc::safemath::safe_acos(std::fabs(evEst.dot(evOldEst))) *
                   (180 / eonc::helpers::pi);
      statsTorque = ewAbsRelErr;
      evOldEst = evEst;
      QUILL_LOG_INFO(log,
                     "[ILanczos] {:9s} {:9s} {:10s} {:14s} {:9.4f} "
                     "{:10.6f} {:7.3f} {:5} n_mobile={}",
                     "----", "----", "----", "----", ew, ewAbsRelErr,
                     statsAngle, i, mobile.size());
      if (ewAbsRelErr < params.lanczos_options().tolerance) {
        QUILL_LOG_INFO(log, "[ILanczos] Tolerance reached: {}",
                       params.lanczos_options().tolerance);
        break;
      }
    } else {
      ew = alpha;
      ewOld = ew;
      evEst = Q.col(0);
      evOldEst = Q.col(0);
      if (lowestEw != 0.0 && params.lanczos_options().quit_early) {
        double Cprev = lowestEw;
        double Cnew = u.dot(Q.col(i));
        ewAbsRelErr = eonc::safemath::safe_div(std::fabs(Cnew - Cprev),
                                               std::fabs(Cprev), 1.0);
        if (ewAbsRelErr <= params.lanczos_options().tolerance) {
          statsAngle = 0.0;
          statsTorque = ewAbsRelErr;
          QUILL_LOG_INFO(log, "[ILanczos] Tolerance reached: {}",
                         params.lanczos_options().tolerance);
          break;
        }
      }
    }

    if (i >= params.lanczos_options().max_iterations - 1) {
      QUILL_LOG_ERROR(log, "[ILanczos] Max iterations");
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

double Lanczos::getEigenvalue() { return lowestEw; }

AtomMatrix Lanczos::getEigenvector() { return lowestEv; }

} // namespace eonc
