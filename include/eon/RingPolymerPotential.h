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

#include "eon/Potential.h"
#include "eon/Tunneling.h"

namespace eonc::tunneling {

/// One structure whose coordinates are the beads of a closed ring. The
/// energy is U_N and the forces are minus its gradient. A call evaluates
/// every bead through the batch potential, so a calculator group can spread
/// the beads. An even ring with the fold set stores one turning point to
/// the other and mirrors the rest.
class RingPolymerPotential : public Potential {
public:
  RingPolymerPotential(BatchPotential beads, long nBeads, long dof,
                       double spring, double energyShift,
                       std::vector<double> eta, bool fold, VectorXd modeDir);

  void force(long nAtoms, const double *positions, const int *atomicNrs,
             double *forces, double *energy, double *variance,
             const double *box) override;

  void forceBatch(long nSystems, long nAtoms, const double *const *positions,
                  const int *const *atomicNrs, double *const *forces,
                  double *energies, double *variances,
                  const double *const *boxes) override;

  [[nodiscard]] bool supportsBatchEvaluation() const noexcept override {
    return true;
  }

  [[nodiscard]] bool requiresIsolatedMoleculeLayout() const noexcept override {
    return true;
  }

  /// The bead batch is a caller-supplied function and is not safe to enter
  /// from two threads.
  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }

  [[nodiscard]] long structureAtoms() const { return nAtoms_; }
  [[nodiscard]] long activeBeads() const { return nActive_; }

  void packActive(const std::vector<VectorXd> &active, double *positions) const;
  [[nodiscard]] std::vector<VectorXd> unpack(const double *positions) const;
  void packMode(const VectorXd &dir, AtomMatrix &mode) const;

private:
  void fixLayout(long nBeads, long dof, bool fold, const VectorXd &modeDir);
  [[nodiscard]] double writeForces(const std::vector<VectorXd> &grad,
                                   const double *positions,
                                   double *forces) const;

  BatchPotential beads_;
  long nBeads_{0};
  long dof_{0};
  long nActive_{0};
  long atomsPerBead_{0};
  int usedOnLast_{3};
  long nAtoms_{0};
  double spring_{0.0};
  double energyShift_{0.0};
  std::vector<double> eta_;
  bool fold_{false};
  VectorXd modeDir_;
};

} // namespace eonc::tunneling
