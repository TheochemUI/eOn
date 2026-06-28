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
// Leng et al. JCP 138, 094110 (2013) Locally Optimal Rotation (Algorithm I).
// FD Hessian-vector products match Lanczos energy-Hessian convention:
//   H v ≈ -(F(x+δv) - F(x)) / δ  with F = forces.
// At most one *new* force evaluation per rotation iteration (F at R0+δ·dir);
// center forces F(R0) are cached once for the whole LOR call.

#include "LORRotation.h"
#include "HelperFunctions.h"
#include "SafeMath.h"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace eonc::helpers;

LORRotation::LORRotation(std::shared_ptr<Matter> matter,
                         const Parameters &params,
                         std::shared_ptr<Potential> pot)
    : LowestEigenmode(pot, params) {
  auto x1Pot = (pot->needsPerImageInstance() && params.main_options.parallel)
                   ? eonc::helpers::makePotential(params)
                   : pot;
  x0 = std::make_shared<Matter>(pot, params);
  x1 = std::make_shared<Matter>(x1Pot, params);
  *x0 = *matter;
  *x1 = *matter;
  totalForceCalls = 0;
  statsRotations = 0;
  eigenvector.resize(3 * matter->numberOfAtoms());
  eigenvector.setZero();
}

VectorXd
LORRotation::translateHUnitOrthoP3(const VectorXd &N, const VectorXd &Theta,
                                   const VectorXd &P, const VectorXd &HN,
                                   const VectorXd &HTheta, const VectorXd &HP) {
  // P_ortho = (I - N Nᵀ - Θ Θᵀ) P  with {N, Θ} orthonormal
  const VectorXd P_ortho = P - N.dot(P) * N - Theta.dot(P) * Theta;
  const double pNrm = P_ortho.norm();
  if (pNrm < 1e-14) {
    return VectorXd::Zero(P.size());
  }
  // H P_ortho = HP - (N·P) HN - (Θ·P) HΘ  (linearity of H — not GS on HP)
  const VectorXd HPortho = HP - N.dot(P) * HN - Theta.dot(P) * HTheta;
  return HPortho / pNrm; // = H (P_ortho / ||P_ortho||) = H P3
}

// Member used only via compute(); keep signature for header compatibility.
VectorXd LORRotation::hessianVector(const VectorXd & /*unused*/,
                                    const VectorXd &x0_r, const VectorXd &v,
                                    const VectorXd &freeMask, double delta) {
  VectorXd F0 = x0->getForcesV();
  VectorXd dir = v.array() * freeMask.array();
  const double nrm = dir.norm();
  if (nrm < 1e-14) {
    return VectorXd::Zero(v.size());
  }
  dir /= nrm;
  x1->setPositionsV(x0_r + delta * dir);
  VectorXd F1 = x1->getForcesV();
  totalForceCalls += 1;
  VectorXd Hv = -(F1 - F0) / delta;
  Hv = Hv.array() * freeMask.array();
  return Hv * nrm;
}

void LORRotation::compute(std::shared_ptr<Matter> matter,
                          AtomMatrix initialDirectionAtomMatrix) {
  const int dim = static_cast<int>(3 * matter->numberOfAtoms());
  const VectorXd freeMask = matter->getFreeV();

  VectorXd N = VectorXd::Map(initialDirectionAtomMatrix.data(), dim);
  N = N.array() * freeMask.array();
  if (N.norm() < 1e-10) {
    N.setRandom();
    N = N.array() * freeMask.array();
  }
  N.normalize();

  *x0 = *matter;
  *x1 = *matter;
  const VectorXd x0_r = x0->getPositionsV();
  const double delta = params.main_options.finiteDifference;

  const int rotmax =
      std::clamp(static_cast<int>(params.dimer_options.rotations_max > 0
                                      ? params.dimer_options.rotations_max
                                      : 20),
                 1, 50);

  const double residualTol = std::max(1e-3, params.dimer_options.torque_min);
  auto relativeResidual = [](double fnorm, double cn) {
    return fnorm / (std::abs(cn) + 1.0);
  };

  curvatureHistory.clear();
  convergedOnResidual = false;
  double bestCN = std::numeric_limits<double>::infinity();
  VectorXd bestN = N;
  VectorXd bestHN = VectorXd::Zero(dim);
  statsRotations = 0;
  totalForceCalls = 0;

  // Cache center forces once (paper: FD products reuse F(R0)).
  const VectorXd F0 = x0->getForcesV();
  totalForceCalls += 1;

  auto applyMask = [&](VectorXd &v) { v = v.array() * freeMask.array(); };

  auto unitize = [&](VectorXd &v) -> double {
    applyMask(v);
    const double n = v.norm();
    if (n > 1e-14) {
      v /= n;
    }
    return n;
  };

  // One new force at R0+δ·dir; F0 cached (Algorithm I: ≤1 new FD per rotation).
  auto hessianAlong = [&](const VectorXd &v) -> VectorXd {
    VectorXd dir = v.array() * freeMask.array();
    const double nrm = dir.norm();
    if (nrm < 1e-14) {
      return VectorXd::Zero(dim);
    }
    dir /= nrm;
    x1->setPositionsV(x0_r + delta * dir);
    const VectorXd F1 = x1->getForcesV();
    totalForceCalls += 1;
    VectorXd Hv = -(F1 - F0) / delta;
    applyMask(Hv);
    return Hv * nrm;
  };

  auto trackBest = [&](double cn, const VectorXd &nVec, const VectorXd &hnVec) {
    if (cn < bestCN) {
      bestCN = cn;
      bestN = nVec;
      bestHN = hnVec;
    }
  };
  // Record Ritz C only when non-increasing (paper quadratic + translation).
  // Strictly non-increasing (1e-4 float). Returns whether the sample was kept
  // (callers must only trackBest on kept samples so ev == min(history)).
  auto appendHistory = [&](double cn) -> bool {
    if (curvatureHistory.empty() || cn <= curvatureHistory.back() + 1e-4) {
      curvatureHistory.push_back(cn);
      return true;
    }
    return false;
  };

  // --- Algorithm I start: H N, F_⊥ ---
  VectorXd HN = hessianAlong(N);
  applyMask(HN);
  double CN = N.dot(HN);
  if (appendHistory(CN)) {
    trackBest(CN, N, HN);
  }

  VectorXd F = HN - CN * N;
  applyMask(F);
  double Fnorm = F.norm();

  QUILL_LOG_INFO(
      log, "[LOR] iter=0 ||F_perp||={:.6e} C_N={:.6f} (Algorithm I start)",
      Fnorm, CN);

  if (relativeResidual(Fnorm, CN) < residualTol) {
    convergedOnResidual = true;
    eigenvalue = CN;
    eigenvector = N;
    QUILL_LOG_INFO(log, "[LOR] converged on residual at start");
    return;
  }

  VectorXd Theta = F / Fnorm;
  VectorXd HTheta = hessianAlong(Theta); // iteration-1 second FD (H Θ)
  applyMask(HTheta);

  // 2×2 Ritz in span{N, Θ}
  Eigen::Matrix2d A2;
  A2(0, 0) = N.dot(HN);
  A2(0, 1) = N.dot(HTheta);
  A2(1, 0) = A2(0, 1);
  A2(1, 1) = Theta.dot(HTheta);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es2(A2);
  Eigen::Vector2d coeffs = es2.eigenvectors().col(0);
  double a = coeffs(0);
  double b = coeffs(1);

  VectorXd Nlin = a * N + b * Theta;
  const double nN = unitize(Nlin);
  VectorXd HNlin = a * HN + b * HTheta;
  if (nN > 1e-14) {
    HNlin /= nN;
  }
  N = Nlin;
  HN = HNlin;

  VectorXd P = Theta;
  VectorXd HP = HTheta;

  CN = N.dot(HN);
  if (appendHistory(CN)) {
    trackBest(CN, N, HN);
  }
  F = HN - CN * N;
  applyMask(F);
  Fnorm = F.norm();
  statsRotations = 1;

  QUILL_LOG_INFO(
      log, "[LOR] iter=1 (2x2) ||F_perp||={:.6e} C_N={:.6f} a={:.4f} b={:.4f}",
      Fnorm, CN, a, b);

  for (int k = 2; k <= rotmax; ++k) {
    if (relativeResidual(Fnorm, CN) < residualTol) {
      convergedOnResidual = true;
      QUILL_LOG_INFO(log,
                     "[LOR] converged residual iter={} ||F||={:.6e} rel={:.6e}",
                     k, Fnorm, relativeResidual(Fnorm, CN));
      break;
    }

    Theta = F / Fnorm;
    // Exactly one new FD this iteration: H · Θ (force translation for N, P)
    HTheta = hessianAlong(Theta);
    applyMask(HTheta);

    // Orthonormalize P vs N,Θ for stable 3×3 (intentional vs paper GEP on B≠I)
    VectorXd P3 = P - N.dot(P) * N - Theta.dot(P) * Theta;
    applyMask(P3);
    const double pNrm = P3.norm();

    auto applyRitz2 = [&]() {
      Eigen::Matrix2d A2b;
      A2b(0, 0) = N.dot(HN);
      A2b(0, 1) = N.dot(HTheta);
      A2b(1, 0) = A2b(0, 1);
      A2b(1, 1) = Theta.dot(HTheta);
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es2b(A2b);
      Eigen::Vector2d c2 = es2b.eigenvectors().col(0);
      a = c2(0);
      b = c2(1);
      VectorXd Nnew = a * N + b * Theta;
      VectorXd HNnew = a * HN + b * HTheta;
      const double nn = unitize(Nnew);
      if (nn > 1e-14) {
        HNnew /= nn;
      }
      const double cNew = Nnew.dot(HNnew);
      if (!curvatureHistory.empty() && cNew > curvatureHistory.back() + 1e-4) {
        QUILL_LOG_INFO(log,
                       "[LOR] curvature stall (2x2 fallback) iter={} "
                       "C_prev={:.6f} C_new={:.6f} (reject)",
                       k, curvatureHistory.back(), cNew);
        N = bestN;
        HN = bestHN;
        CN = bestCN;
        F = HN - CN * N;
        applyMask(F);
        Fnorm = F.norm();
        statsRotations = k;
        return false;
      }
      P = Theta;
      HP = HTheta;
      N = Nnew;
      HN = HNnew;
      CN = cNew;
      if (appendHistory(CN)) {
        trackBest(CN, N, HN);
      }
      F = HN - CN * N;
      applyMask(F);
      Fnorm = F.norm();
      statsRotations = k;
      return true;
    };

    if (pNrm < 1e-8) {
      if (!applyRitz2()) {
        break;
      }
      continue;
    }
    P3 /= pNrm;
    // H·P3 by force-translation linearity (not ambient GS on HP):
    // H P3 = (HP - (N·P) HN - (Θ·P) HΘ) / ||P_ortho|| with P3 = P_ortho/||...||
    VectorXd HP3 = translateHUnitOrthoP3(N, Theta, P, HN, HTheta, HP);
    applyMask(HP3);
    if (HP3.norm() < 1e-14) {
      if (!applyRitz2()) {
        break;
      }
      continue;
    }

    Eigen::Matrix3d A3 = Eigen::Matrix3d::Zero();
    const VectorXd basis[3] = {N, Theta, P3};
    const VectorXd Hbasis[3] = {HN, HTheta, HP3};
    for (int i = 0; i < 3; ++i) {
      for (int j = i; j < 3; ++j) {
        A3(i, j) = basis[i].dot(Hbasis[j]);
        A3(j, i) = A3(i, j);
      }
    }
    A3 = 0.5 * (A3 + A3.transpose());
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> ges(A3);
    if (ges.info() != Eigen::Success) {
      QUILL_LOG_WARNING(log, "[LOR] 3x3 eigen failed at iter={}; stop", k);
      break;
    }
    const Eigen::Vector3d c3 = ges.eigenvectors().col(0);
    a = c3(0);
    b = c3(1);
    const double c = c3(2);

    // Save pre-update for stall revert without FD
    const VectorXd Nprev = N;
    const VectorXd HNprev = HN;
    const double CNprev = CN;

    VectorXd Nnew = a * N + b * Theta + c * P3;
    VectorXd Pnew = b * Theta + c * P3;
    VectorXd HNnew = a * HN + b * HTheta + c * HP3;
    VectorXd HPnew = b * HTheta + c * HP3;

    const double nNew = unitize(Nnew);
    if (nNew > 1e-14) {
      HNnew /= nNew;
    }
    const double pNew = unitize(Pnew);
    if (pNew > 1e-14) {
      HPnew /= pNew;
    }

    const double CNnew = Nnew.dot(HNnew);
    if (!curvatureHistory.empty() && CNnew > curvatureHistory.back() + 1e-4) {
      QUILL_LOG_INFO(log,
                     "[LOR] curvature stall iter={} C_prev={:.6f} C_new={:.6f} "
                     "(reject update; keep prior mode, no extra FD)",
                     k, curvatureHistory.back(), CNnew);
      N = bestN;
      HN = bestHN;
      CN = bestCN;
      F = HN - CN * N;
      applyMask(F);
      Fnorm = F.norm();
      statsRotations = k;
      break;
    }

    N = Nnew;
    P = Pnew;
    HN = HNnew;
    HP = HPnew;
    CN = CNnew;
    if (appendHistory(CN)) {
      trackBest(CN, N, HN);
    }
    (void)Nprev;
    (void)HNprev;
    (void)CNprev;

    F = HN - CN * N;
    applyMask(F);
    Fnorm = F.norm();
    statsRotations = k;

    QUILL_LOG_INFO(log,
                   "[LOR] iter={} (3x3) ||F_perp||={:.6e} C_N={:.6f} a={:.4f} "
                   "b={:.4f} c={:.4f}",
                   k, Fnorm, CN, a, b, c);

    if (relativeResidual(Fnorm, CN) < residualTol) {
      convergedOnResidual = true;
      break;
    }
  }

  // Return best softest mode from accepted Ritz steps (no mandatory final FD).
  N = bestN;
  HN = bestHN;
  CN = bestCN;
  eigenvalue = CN;
  eigenvector = N;
  QUILL_LOG_INFO(log,
                 "[LOR] done rotations={} force_calls={} C_N={:.6f} "
                 "converged_residual={}",
                 statsRotations, totalForceCalls, eigenvalue,
                 convergedOnResidual);
}
