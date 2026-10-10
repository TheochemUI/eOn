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

// Path-integral quantum transition-state theory (Voth, Chandler and
// Miller, J. Chem. Phys. 91, 7749 (1989)): the centroid potential of mean
// force along a linear mass-weighted coordinate s from ring polymers whose
// centroid is held on parallel planes, and the rate from its value on the
// dividing plane, optionally times the ring-polymer MD transmission factor
// on that plane (Craig and Manolopoulos, J. Chem. Phys. 122, 084106 (2005);
// Suleimanov, Allen and Green, Comput. Phys. Commun. 184, 833 (2013)).

#include "eon/Eigen.h"
#include "eon/ParametersOptions.h"
#include "eon/PathIntegral.h"
#include "eon/Potential.h"

#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace eonc {
class Matter;
class Parameters;
} // namespace eonc

namespace eonc::piqtst {

/// The coordinate s = n . M^(1/2) (x - reference), n a unit vector in
/// mass-weighted coordinates (zero on fixed coordinates). The plane at s is
/// the set of centroids with that value.
struct Coordinate {
  long atoms{0};
  std::vector<double> masses; // per atom, amu
  std::vector<int> numbers;   // per atom
  std::vector<char> free;     // per coordinate, 3 * atoms
  VectorXd reference;         // Cartesian, 3 * atoms, Angstrom
  VectorXd direction;         // mass-weighted unit vector, 3 * atoms
  const double *box{nullptr}; // 3 x 3 cell, or null
};

struct ScanOptions {
  /// Plane positions in amu^0.5 Angstrom, ascending. The last plane is
  /// the dividing surface s*.
  std::vector<double> planes;
  long equilibration{500};
  long production{2000};
  /// Equal blocks of the production run for the standard error.
  long blocks{10};
  /// Beads, temperature, units (kB in eV / K, hbar in eV time units),
  /// time step and thermostat of the ring.
  pathintegral::Options ring;
  /// Cartesian centroid to start the ring at on the plane at s. Empty
  /// starts on the line reference + s M^(-1/2) n.
  std::function<VectorXd(double)> seed;
};

struct Plane {
  double s{0.0};
  /// dF/ds = -<a . f_centroid> / |a|^2, a = M^(1/2) n the plane's
  /// Cartesian normal, eV / (amu^0.5 Angstrom), and its block standard
  /// error. Any u with a . u = 1, M^(-1/2) n among them, has the same mean,
  /// since the force within the plane averages to zero; the normal
  /// component is the estimator.
  double meanForce{0.0};
  double meanForceError{0.0};
  /// F(s) - F(s_0) by the trapezoid rule over the planes, eV, and its
  /// standard error from the mean-force errors.
  double freeEnergy{0.0};
  double freeEnergyError{0.0};
  /// Production average of the centroid, Cartesian.
  VectorXd centroid;
  /// Root-mean-square bead displacement from the centroid per Cartesian
  /// coordinate, Angstrom, averaged over the production run.
  std::vector<double> spread;
  long batches{0};
};

/// Samples one ring per plane and integrates the mean force. One ring is
/// carried from plane to plane: its centroid moves to the next seed and the
/// internal modes keep their state. A step is one force batch over the
/// beads; moving to a plane costs one more, so a plane takes
/// equilibration + production + 1.
std::vector<Plane> scan(Potential &pot, const Coordinate &coordinate,
                        const ScanOptions &options);

/// Trapezoid integral of the mean forces, F(s_0) = 0, with errors from
/// independent planes. Fills freeEnergy and freeEnergyError.
void integrate(std::vector<Plane> &planes);

struct Rate {
  /// Index of the plane with the lowest F before s*, the reactant.
  long reactant{0};
  /// F(s*) - F(reactant), eV, and its error; not positive when s* is no
  /// barrier on the centroid free energy.
  double barrier{0.0};
  double barrierError{0.0};
  /// ln k with k in inverse eOn time units (sqrt(amu Angstrom^2 / eV)),
  /// and its standard error.
  double logRate{0.0};
  double logRateError{0.0};
  /// beta (F(s_0) - F(reactant)): the reactant integral is cut at s_0, so
  /// a small value means the first plane is too close to the well.
  double firstPlaneHeight{0.0};
};

/// k = (1/2) sqrt(2 / (pi beta)) exp(-beta F(s*)) /
///     int_{s_0}^{s*} exp(-beta F(s)) ds,
/// the mass-weighted coordinate carrying unit mass, beta in 1 / eV. The
/// integral is the trapezoid rule over the planes; the last is s*.
Rate rate(const std::vector<Plane> &planes, double beta);

struct RecrossingOptions {
  /// The dividing plane s*, amu^0.5 Angstrom.
  double s{0.0};
  /// Thermostatted steps on the plane before the first parent.
  long equilibration{500};
  /// Parent configurations, each this many thermostatted steps after the
  /// last.
  long parents{0};
  long spacing{50};
  /// Momentum draws per parent; each runs forward and reversed.
  long children{20};
  /// Unconstrained, thermostat-free steps per child of ring.dt.
  long steps{0};
  /// The parents' ring and thermostat. The children use its beads,
  /// temperature and time step.
  pathintegral::Options ring;
  /// Cartesian centroid to start the parent ring at. Empty starts on
  /// reference + s M^(-1/2) n.
  std::function<VectorXd(double)> seed;
};

struct Recrossing {
  /// t = step * ring.dt, from 0 to steps * ring.dt, and
  /// kappa(t) = <sdot(0) h(s(t) - s*)> / <sdot(0) h(sdot(0))>.
  /// At t = 0 the side is that of sdot(0), the t -> 0+ limit, so
  /// kappa(0) = 1.
  std::vector<double> time;
  std::vector<double> kappa;
  /// Mean of kappa(t) over the last quarter of the times, and its
  /// jackknife standard error over parents.
  double plateau{0.0};
  double plateauError{0.0};
  long trajectories{0};
  long batches{0};
};

/// Bennett-Chandler transmission at s*: parents sampled with the centroid
/// held on the plane under PILE, whatever ring.thermostat says, children
/// with Maxwell-Boltzmann ring momenta at beta / P and the plane released,
/// propagated by RPMD. The centroid velocity along s is
/// sdot = n . M^(1/2) v_centroid. Each child run costs steps + 1 force
/// batches: the first force at the parent's beads, then one per step.
Recrossing recrossing(Potential &pot, const Coordinate &coordinate,
                      const RecrossingOptions &options);

/// Throws std::invalid_argument on an inconsistent [Instanton] pi_* key.
void validateOptions(const instanton_options_t &o);

/// [Instanton] mode rate with pi_planes > 0: the planes and rate at each
/// temperature, and with pi_recrossing_parents > 0 the transmission factor
/// on the top plane and k_RPMD = kappa k_PI-QTST. hSaddle and the band in pathQ
/// are in the instanton's mass-weighted coordinates over the free atoms,
/// measured from the reactant; saddle is aligned to the reactant. Appends
/// results.dat keys for the last (lowest) temperature to extras and returns the
/// files written.
std::vector<std::string>
runAfterInstanton(const Parameters &params, Potential &pot,
                  const Matter &reactant, const Matter &saddle,
                  const MatrixXd &hSaddle, const std::vector<VectorXd> &pathQ,
                  const std::vector<double> &temperatures,
                  std::vector<std::pair<std::string, double>> &extras);

} // namespace eonc::piqtst
