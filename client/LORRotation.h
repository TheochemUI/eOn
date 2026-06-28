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
#pragma once
#include "Eigen.h"
#include "EonLogger.h"
#include "LowestEigenmode.h"
#include "Matter.h"
#include "Parameters.h"

namespace eonc {

// Leng et al., J. Chem. Phys. 138, 094110 (2013) — Locally Optimal Rotation
// (LOR) Algorithm I for softest-mode finding in the dimer rotation step. At
// most one new force evaluation per rotation iteration via FD H·v + force
// translation of prior H·N / H·F_⊥ / H·P products.
class LORRotation : public LowestEigenmode {
private:
  eonc::log::Scoped log;
  VectorXd eigenvector;
  double eigenvalue{0.0};
  std::shared_ptr<Matter> x0;
  std::shared_ptr<Matter> x1;

  // FD Hessian-vector product at x0 along unit direction v (free atoms only).
  // Uses gradients g = -forces; H v ≈ (g(x0+δv) - g(x0)) / δ.
  VectorXd hessianVector(const VectorXd &g0, const VectorXd &x0_r,
                         const VectorXd &v, const VectorXd &freeMask,
                         double delta);

public:
  LORRotation(std::shared_ptr<Matter> matter, const Parameters &params,
              std::shared_ptr<Potential> pot);
  ~LORRotation() = default;

  void compute(std::shared_ptr<Matter> matter, AtomMatrix initialDirection);
  double getEigenvalue() { return eigenvalue; }
  AtomMatrix getEigenvector() {
    return AtomMatrix::Map(eigenvector.data(), eigenvector.size() / 3, 3);
  }

  // Last LOR curvature sequence (for unit tests / monotonicity checks)
  std::vector<double> curvatureHistory;
  bool convergedOnResidual{false};

  /// Force-translation identity for H·P3 when P3 is the *unit* Gram-Schmidt
  /// residual of P against orthonormal {N, Θ} and H is linear:
  ///   P_ortho = P - (N·P) N - (Θ·P) Θ
  ///   P3 = P_ortho / ||P_ortho||
  ///   H P3 = (H P - (N·P) H N - (Θ·P) H Θ) / ||P_ortho||
  /// (Not Gram-Schmidt on HP in ambient space.) Returns zero vector if
  /// ||P_ortho|| is negligible. Shipped LOR 3×3 path uses this helper.
  [[nodiscard]] static VectorXd
  translateHUnitOrthoP3(const VectorXd &N, const VectorXd &Theta,
                        const VectorXd &P, const VectorXd &HN,
                        const VectorXd &HTheta, const VectorXd &HP);
};

} // namespace eonc

using eonc::LORRotation;
