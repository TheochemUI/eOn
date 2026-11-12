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
#include "ObjectiveFunction.h"
#include "Parameters.h"

#include <array>
#include <iosfwd>
#include <memory>
#include <string_view>
#include <vector>

namespace eonc {
class JobResultEnvelope;
}

namespace eonc::solid {

/// Lower-triangular strain components, in order xx, yx, yy, zx, zy, zz.
inline constexpr int kStrainRows[6] = {0, 1, 1, 2, 2, 2};
inline constexpr int kStrainCols[6] = {0, 0, 1, 0, 1, 2};

struct Options {
  bool enabled{false};
  double weight{1.0};
  /// Cauchy stress added to the potential stress before the cell force,
  /// in eV/Angstrom^3. Hydrostatic pressure is already on the diagonal.
  Matrix3d external{Matrix3d::Zero()};
  /// True when that strain component is a degree of freedom.
  std::array<bool, 6> mask{{true, true, true, true, true, true}};
  /// Largest |Δε| of one component in a single coordinate update.
  double maxStrainStep{0.05};
  /// Largest |ε| relative to the starting cell. Past this the update throws.
  double maxStrain{0.20};
};

/// How a converged solid-state dimer checks that the stationary point
/// is index 1. The defaults match the check used on the copper shear
/// rows: eight lowest modes, a central step of 1e-4, and a negative
/// eigenvalue below -1e-4. Six escapes of length 0.02 leave a
/// higher-index point along its softest negative mode.
struct IndexCheck {
  int lowest{8};
  double negative{-1.0e-4};
  double step{1.0e-4};
  int escapes{6};
  double escape{0.02};
  /// Rigid translations of the free atoms are removed before the count.
  bool project{true};
};

/// Lowest projected eigenvalues, ascending. Eigenvectors are columns.
struct IndexSpectrum {
  int negative{0};
  std::vector<double> values;
  MatrixXd vectors;
  /// Generalized evaluations spent on the Hessian actions.
  int evaluations{0};

  [[nodiscard]] bool isFirstOrder() const { return negative == 1; }
};

[[nodiscard]] Options optionsFrom(const Parameters &params);

/// True when the cell has at least one free strain component.
/// The frozen cell (enabled off, or a mask of all zeros) is not this.
[[nodiscard]] bool movesCell(const Parameters &params);

/// Reference chart for one minimization or saddle search.
///
/// The cell is h = h0 (I + ε) with ε lower triangular. The atomic
/// coordinates are the Cartesian positions in the reference cell, so a
/// pure cell step keeps fractional coordinates fixed. The scaled strain
/// is q = J ε with J = V0^(1/3) N^(1/6) weight, held at the reference
/// volume. doi:10.1063/1.3684549
class Chart {
public:
  Chart(Matter &matter, const Options &options, AtomMatrix *mode = nullptr);

  [[nodiscard]] int degreesOfFreedom() const {
    return 3 * static_cast<int>(free_.size()) +
           static_cast<int>(active_.size());
  }
  [[nodiscard]] double jacobian() const { return jacobian_; }
  [[nodiscard]] const Matrix3d &referenceCell() const { return h0_; }
  [[nodiscard]] const std::vector<long> &freeAtoms() const { return free_; }
  [[nodiscard]] const std::vector<int> &activeStrain() const { return active_; }
  /// >1 shrinks a stiff cell direction in the joint chart.
  [[nodiscard]] double cellScale(size_t index) const {
    return cellScale_.at(index);
  }

  [[nodiscard]] VectorXd coordinates(const Matter &matter) const;
  void setCoordinates(Matter &matter, const VectorXd &x) const;

  /// Scale stiff cell directions so their curvature matches an atomic one.
  /// The Jacobian metric is unchanged: this only rescales the cell block
  /// the min-mode rotation and the optimizer step in.
  void preconditionCell(Matter &matter);

  /// dE/dq in the chart. One force evaluation when the geometry is dirty.
  [[nodiscard]] VectorXd energyGradient(Matter &matter) const;

  [[nodiscard]] VectorXd generalizedForce(Matter &matter) const {
    return -energyGradient(matter);
  }

  /// Norm of the force in the Jacobian metric, before the cell-block scale.
  /// Convergence uses this, so a stiff direction cannot hide a real force.
  [[nodiscard]] double physicalResidual(Matter &matter) const;

  /// `metric` is the optimizer convergence_metric. The cell block is
  /// measured before the diagonal scale.
  [[nodiscard]] double convergence(Matter &matter,
                                   const std::string &metric) const;

  /// Minimum-image difference of two coordinate vectors. The cell block
  /// is not periodic.
  [[nodiscard]] VectorXd difference(const VectorXd &a, const VectorXd &b) const;

  /// Atomic rows of a joint vector. Fixed atoms are zero.
  [[nodiscard]] AtomMatrix unpackAtoms(const VectorXd &joint) const;
  /// Lower-triangular strain direction, not scaled by J.
  [[nodiscard]] Matrix3d unpackStrain(const VectorXd &joint) const;
  [[nodiscard]] VectorXd packAtoms(const AtomMatrix &mode) const;

  /// Two smallest curvatures in the joint space. A first-order saddle has
  /// a negative lowest value and a positive second one.
  [[nodiscard]] std::array<double, 2>
  lowestCurvatures(Matter &matter, VectorXd direction,
                   const Parameters &params) const;

  /// Diagonal scale on the cell block. The index check clears it so the
  /// Hessian is the physical Jacobian metric.
  [[nodiscard]] std::vector<double> cellScales() const { return cellScale_; }
  void setCellScales(std::vector<double> scales);

  /// E + V0 σ_ext : ε. The upper triangle of ε does no work.
  [[nodiscard]] double enthalpy(Matter &matter) const;

private:
  Options options_;
  Matrix3d h0_{Matrix3d::Identity()};
  AtomMatrix r0_;
  double jacobian_{1.0};
  double referenceVolume_{1.0};
  std::vector<long> free_;
  std::vector<int> active_;
  /// Multiplier on q = J ε. Stiff components are > 1.
  std::vector<double> cellScale_;

  [[nodiscard]] Matrix3d deformation(const Matter &matter) const;
  void readStrain(const Matrix3d &deformation, Matrix3d &eps) const;
};

[[nodiscard]] std::shared_ptr<ObjectiveFunction>
makeVariableCellObjective(Matter &matter, const Parameters &params);

struct ModeEstimate {
  VectorXd mode;
  double curvature{0.0};
};

[[nodiscard]] ModeEstimate
lowestMode(Matter &matter, const Chart &chart, VectorXd direction,
           const Parameters &params);

/// One finished solid-state min-mode search. The chart stays valid so a
/// caller can step off the saddle along the joint mode.
struct SearchState {
  Chart chart;
  VectorXd mode;
  double curvature{0.0};
  int status{0};
  int iteration{0};
  int forceCalls{0};
  /// Potential calls spent building projected Hessians and walking off
  /// a higher-index point. Already included in `forceCalls`.
  int indexForceCalls{0};
  /// Negative projected eigenvalues at an accepted or rejected point.
  /// -1 when the check did not run.
  int negativeModes{-1};
  /// The force tolerance was met and the projected Hessian is not index 1.
  bool rejected{false};

  explicit SearchState(Chart chartIn) : chart(std::move(chartIn)) {}
};

/// Projected spectrum of the joint Hessian at the current geometry.
/// The matter is put back when the differences are done.
[[nodiscard]] IndexSpectrum indexAt(Matter &matter, const Chart &chart,
                                    const IndexCheck &check = {});

/// Same count on a frozen cell: the atomic Hessian, with rigid
/// translations projected out when `check.project` is set.
[[nodiscard]] IndexSpectrum indexAtoms(Matter &matter,
                                       const IndexCheck &check = {});

/// Step off a higher-index point along its softest negative mode, on
/// the downhill side, until that curvature is no longer below the
/// cutoff. `hardMode` is the lowest eigenvector. False when neither
/// sign descends, or the walk returns to a point already left.
[[nodiscard]] bool leaveHigherIndex(Matter &matter, const Chart &chart,
                                    VectorXd &hardMode,
                                    const IndexSpectrum &spectrum,
                                    const IndexCheck &check,
                                    std::vector<VectorXd> &rejectedAt);

[[nodiscard]] bool
leaveHigherIndexAtoms(Matter &matter, VectorXd &hardMode,
                      const IndexSpectrum &spectrum, const IndexCheck &check,
                      std::vector<VectorXd> &rejectedAt);

[[nodiscard]] std::unique_ptr<SearchState>
runMinMode(Matter &matter, AtomMatrix mode, double reactantEnergy,
           const Parameters &params, long maxIterations);

/// `image` starts at the saddle. `scale` is a step in the joint chart.
void placeAlongMode(const SearchState &state, Matter &image, double scale);

/// Cauchy stress, cell volume, and the elastic dipole P = V σ.
/// Nothing is written when the potential does not report a stress.
void appendElasticDipole(JobResultEnvelope &env, Matter &matter,
                         std::string_view prefix);
void writeElasticDipole(std::ostream &out, Matter &matter,
                        std::string_view prefix);

} // namespace eonc::solid
