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

#include "FiniteDifference.h"
#include "Matter.h"
#include "Parameters.h"

#include <vector>

namespace eonc {

/// Whether `removed` zero-frequency modes are what the structure's
/// symmetries give: 6 for a free cluster (3 translations and 3 rotations), 5
/// for a linear one, 3 for a periodic cell, where rotations are no symmetry
/// of the lattice, and none once any atom is fixed.
bool trivialModeCountIsPhysical(long removed, long fixedAtoms);

/// 1 eV in cm^-1, e / (h c) from the exact SI values.
inline constexpr double kEvToWavenumber = 8065.543937349212;

/// Cartesian displacement of every atom along one mass-weighted mode over
/// the mobile degrees of freedom `atoms`: x = q / sqrt(m), unit norm, zero
/// on atoms outside `atoms`. Row-major N x 3.
std::vector<double> cartesianMode(const Matter &matter, const VectorXi &atoms,
                                  const Eigen::Ref<const VectorXd> &mode);

/// Writes one frame of `matter` per mode to `path`, the mode as the
/// displacements section and `mode_eigenvalue` (eV / (Angstrom^2 amu)),
/// `hbar_omega` (eV, negative for an imaginary mode) and `wavenumber`
/// (cm^-1, same sign) as frame metadata.
bool writeNormalModes(Matter &matter, const VectorXi &atoms,
                      const VectorXd &eigenvalues, const MatrixXd &modes,
                      const std::string &path);

class Hessian {
public:
  Hessian(const Parameters &params, Matter *matter);
  ~Hessian() = default;

  MatrixXd getHessian(Matter *matterIn, const VectorXi &atomsIn);
  VectorXd getFreqs(Matter *matterIn, const VectorXi &atomsIn);
  VectorXd removeZeroFreqs(const VectorXd &freqs);
  /// Eigenvectors of the mass-weighted Hessian, one column per eigenvalue
  /// of getFreqs(), over the mobile degrees of freedom. Empty unless
  /// [Hessian] write_modes is set.
  [[nodiscard]] const MatrixXd &getModes() const noexcept { return modes; }
  /// Whether a finished Hessian goes to hessian.dat (on by default). A
  /// caller that takes many Hessians, one per instanton bead, turns it off.
  void writeHessianFile(bool on) noexcept { writeFile = on; }

private:
  Matter *matter;
  const Parameters &parameters;

  MatrixXd hessian;
  VectorXd freqs;
  MatrixXd modes;
  bool writeFile = true;

  VectorXi atoms;
  bool calculate();
  bool finalizeHessian(int size);
  bool calculateColored(double cutoff, double dr, FdScheme scheme);
  bool calculateSerial(double dr, FdScheme scheme);
  bool calculateBatched(double dr, FdScheme scheme);
  eonc::log::Scoped log;
};

/// Greedy coloring of an undirected graph. ``adj[v]`` lists neighbors of
/// ``v`` (self-loops ignored). Colors are dense from 0 in vertex order.
std::vector<int>
greedyColorCutoffGraph(const std::vector<std::vector<int>> &adj);

/// Conflict graph on ``atoms``: two mobile atoms share an edge when their
/// closed cutoff neighborhoods intersect, then greedy-colored. Empty when
/// the neighbor list cannot be built.
std::vector<int> colorMobileCutoffGraph(const Matter &matter,
                                        const VectorXi &atoms, double cutoff);

} // namespace eonc
