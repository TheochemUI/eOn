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
#include <functional>
#include <string>

#include "eon/Eigen.h"

namespace eonc {

/**
 * \brief Colored-noise (generalized Langevin) thermostat.
 *
 * Ceriotti, Bussi & Parrinello, "Colored-noise thermostats a la
 * carte", J. Chem. Theory Comput. 6, 1170 (2010): each Cartesian
 * degree of freedom carries `n_aux` auxiliary momenta, and the
 * extended momentum vector evolves under the exact Ornstein-Uhlenbeck
 * map
 *
 *   (p, s) <- T (p, s) + S xi,   T = exp(-A dt/2),
 *   S S^T = C - T C T^T,         C = kB T * I (canonical sampling),
 *
 * applied as half-steps around a velocity Verlet integrator. The
 * frequency-dependent friction shaped by A thermalises every mode of
 * a stiff spectrum at comparable efficiency, where white noise
 * critically damps only one band.
 *
 * The drift matrix A is an operator input ((n_aux+1) x (n_aux+1)
 * plain-text rows, inverse internal time units, the layout produced
 * by the gle4md generator). A 1x1 matrix [gamma] degenerates to
 * standard white-noise Langevin, so the thermostat carries no hidden
 * per-system tuning of its own.
 *
 * All degrees of freedom propagate together: with the extended state
 * held as a (1 + n_aux) x n_dof matrix Z (row 0 the mass-scaled
 * velocities, rows 1.. the auxiliaries), one half-step is
 * Z <- T Z + S Xi with Xi standard normal.
 */
class GleThermostat {
public:
  /// Build from the drift matrix, target kB*T (eV), the half-step dt
  /// (internal units), and the number of free degrees of freedom.
  GleThermostat(const MatrixXd &a_drift, double kbt, double dt_half,
                long n_dof);

  /// Load a drift matrix from a gle4md-layout text file ('#' starts a
  /// comment). Returns an empty (0 x 0) matrix on parse failure.
  static MatrixXd loadDriftMatrix(const std::string &path);

  /// One OU half-step on the mass-scaled velocities; `gauss` supplies
  /// standard-normal draws from the job's deterministic stream.
  void apply(VectorXd &vel, const VectorXd &masses3N,
             const std::function<double()> &gauss);

  /// True when construction produced a usable propagator.
  bool valid() const { return m_valid; }

private:
  bool m_valid{false};
  MatrixXd m_T; //!< exp(-A dt/2)
  MatrixXd m_S; //!< noise factor, S S^T = kbt (I - T T^T)
  MatrixXd m_Z; //!< extended state, (1 + n_aux) x n_dof
};

} // namespace eonc

using eonc::GleThermostat;
