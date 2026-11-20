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
#include "Matter.h"
#include "Parameters.h"
#include "RingDeterminant.h"
#include "Tunneling.h"

#include <array>
#include <atomic>
#include <memory>
#include <stdexcept>
#include <vector>

namespace eonc {

/// Mass-weighted coordinates over the free atoms, measured from a reference
/// structure under its minimum image.
class MassWeighted {
public:
  explicit MassWeighted(const Matter &reference);

  long dimension() const;
  VectorXi freeAtoms() const;
  std::vector<double> spreadAbout(const std::vector<VectorXd> &beads,
                                  VectorXd &centroid, long atoms) const;
  VectorXd toQ(const Matter &matter) const;
  void place(const VectorXd &q, Matter &matter) const;
  MatrixXd rigidGenerators(const Matter &matter) const;
  void markRotationZeroModes(const MatrixXd &hessian,
                             const MatrixXd &generators,
                             std::array<bool, 3> &keep,
                             std::array<double, 3> &residual) const;
  MatrixXd rigidBasis(const Matter &matter,
                      const std::array<bool, 3> &rotations) const;
  const std::vector<double> &sqrtMasses() const;
  double minimumSqrtMass() const;
  VectorXd referenceFree() const;
  VectorXd gradient(const AtomMatrix &forces) const;

private:
  const Matter &reference_;
  std::vector<long> free_;
  std::vector<double> sqrtMass_;
};

/// Remove the mass-weighted translation and, for a finite free system, the
/// best proper rotation from Kabsch RMSD alignment.
void alignRigid(const Matter &reference, Matter &matter);

/// The potential has no exact mass-weighted Hessian-vector product where a
/// caller explicitly requires one.
class HessianVectorProductUnavailableError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

/// Shared potential, coordinate, Hessian-vector product, and finite-difference
/// Hessian adapter used by tunnelling jobs.
class TunnelingJobContext {
public:
  struct Options {
    bool requireWholeAtomConstraints{false};
    bool forwardNonPeriodicCell{false};
    /// Independent potential instances used by each ring-Hessian product.
    /// One preserves the serial job and benchmark path.
    long hessianProductThreads{1};
  };

  TunnelingJobContext(const Parameters &parameters,
                      std::shared_ptr<Potential> potential,
                      const Matter &reference);
  TunnelingJobContext(const Parameters &parameters,
                      std::shared_ptr<Potential> potential,
                      const Matter &reference, Options options);

  const MassWeighted &coordinates() const;
  void evaluate(const std::vector<VectorXd> &q, std::vector<double> &energies,
                std::vector<VectorXd> &gradients);
  MatrixXd hessian(const VectorXd &q, bool projectRotations);
  MatrixXd hessianAt(const VectorXd &q);
  MatrixXd beadHessianAt(const VectorXd &q);
  struct LowestMode {
    double curvature{0.0};
    VectorXd vector;
    double residual{0.0};
    long iterations{0};
  };
  /// The lowest mass-weighted Hessian mode at q by Lanczos, using an exact
  /// product when available and the declared finite-difference fallback
  /// otherwise, to the Lanczos tolerance relative to the curvature.
  LowestMode lowestModeAt(const VectorXd &q, const VectorXd &seed);
  /// H(q) v from the potential's exact product when available, or from one
  /// finite-difference Hessian assembled and reported at q otherwise.
  tunneling::HessianProduct hessianProductAt(const VectorXd &q);
  tunneling::RingHessianProducts
  ringHessianProducts(const std::vector<VectorXd> &beads);
  [[nodiscard]] bool
  usesFiniteDifferenceHessianFallback() const noexcept;
  /// Products taken through the closures above and lowestModeAt.
  long hessianVectorProducts() const;
  /// Marks which rotations are zero modes at q, by the rule of
  /// tunneling::rotationZeroModes, from Rayleigh quotients and a Lanczos
  /// softest vibration on the complement of the rigid generators.
  void findRotationZeroModes(const VectorXd &q);
  const std::array<bool, 3> &rotationZero() const;
  const std::array<double, 3> &rotationResidual() const;

private:
  void reportFiniteDifferenceFallback();
  LowestMode lowestMode(const tunneling::HessianProduct &product,
                        const VectorXd &seed, const MatrixXd &excluded);

  const Parameters &parameters_;
  std::shared_ptr<Potential> potential_;
  Matter reference_;
  MassWeighted coordinates_;
  std::vector<std::unique_ptr<Matter>> pool_;
  std::array<bool, 3> rotationZero_{{false, false, false}};
  std::array<double, 3> rotationResidual_{{0.0, 0.0, 0.0}};
  bool rotationsKnown_{false};
  bool reportedFiniteDifferenceFallback_{false};
  bool forwardNonPeriodicCell_{false};
  long hessianProductThreads_{1};
  std::shared_ptr<std::atomic<long>> products_{
      std::make_shared<std::atomic<long>>(0)};
};

} // namespace eonc
