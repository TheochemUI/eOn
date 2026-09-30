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

// One-dimensional tunnelling along a mass-weighted minimum energy path.
//
// A two-level system in a glass is a pair of adjacent minima that the system
// tunnels between at about one kelvin. Its splitting follows from a WKB
// integral along the mass-weighted path joining the minima (Damart and
// Rodney; Khomenko et al.; Mocanu et al.). Energies are in eV, lengths in
// Angstrom and masses in amu, so mass-weighted lengths are in
// amu^0.5 Angstrom and hbar is kHbar.

#include "Matter.h"

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace eonc::tunneling {

/// hbar in eV^0.5 amu^0.5 Angstrom, correctly rounded from the exact SI
/// values (h = 6.62607015e-34 J s, e = 1.602176634e-19 C) and the CODATA
/// 2022 dalton, 1.66053906892e-27 kg.
inline constexpr double kHbar = 0.06465415129579072;

/// Boltzmann constant in eV / K, correctly rounded from the exact
/// 1.380649e-23 J / K.
inline constexpr double kBoltzmann = 8.617333262145177e-5;

/// sqrt(sum_i m_i |b_i - a_i|^2) under the minimum image of a's cell.
double massWeightedDistance(const Matter &a, const Matter &b);

/// Cumulative mass-weighted arc length at each image of a band.
std::vector<double>
massWeightedPath(const std::vector<std::shared_ptr<Matter>> &band);

/// Energy along the path, interpolated as a monotone cubic (Fritsch and
/// Carlson), flat at both ends because the ends of a band are minima. A
/// monotone interpolant cannot invent a dip below the data, which would
/// change the forbidden region the WKB integral runs over.
class Profile {
public:
  Profile(std::vector<double> s, std::vector<double> v);
  double operator()(double x) const;
  const std::vector<double> &s() const { return s_; }
  const std::vector<double> &v() const { return v_; }

private:
  std::vector<double> s_, v_, m_;
};

/// d2V/ds2 at one end of the path, in eV / (amu Angstrom^2), from a least
/// squares fit of a s^2 + b s^3 to the points within half the barrier above
/// that end. The cubic term takes the anharmonic part a parabola alone would
/// fold into the curvature.
double wellCurvature(const Profile &p, bool leftEnd);

/// hbar omega in eV for a mass-weighted curvature.
double hbarOmega(double curvature);

/// (1/hbar) integral sqrt(2 (V(s) - E)) ds over the path where V > E.
double wkbAction(const Profile &p, double energy, int points = 4001);

struct Splitting {
  double delta = 0.0;           ///< product minus reactant, eV
  double barrier = 0.0;         ///< path maximum above the reactant, eV
  double hwReactant = 0.0;      ///< hbar omega of the reactant well, eV
  double hwProduct = 0.0;       ///< hbar omega of the product well, eV
  double referenceEnergy = 0.0; ///< the level tunnelling happens at, eV
  double action = 0.0;          ///< the WKB exponent, dimensionless
  double delta0 = 0.0;          ///< tunnelling splitting, eV
  /// Both barriers stand above hbar omega; below that WKB is not the right
  /// tool and the number is reported but flagged.
  bool deepWells = false;
  double tlsEnergy() const; ///< sqrt(delta^2 + delta0^2), eV
};

/// delta0 = (hbar omega / pi) exp(-S) with the Landau and Lifshitz
/// prefactor. The level is the higher of the two harmonic ground states,
/// V_i + hbar omega_i / 2, and omega is the geometric mean of the wells;
/// for a symmetric double well both reduce to the textbook formula.
Splitting wkbSplitting(const Profile &p, double hwReactant, double hwProduct);

/// The splitting of a converged band, with the well frequencies from the
/// band's curvature at each end.
Splitting bandSplitting(const std::vector<std::shared_ptr<Matter>> &band,
                        double referenceEnergy);

// Ring-polymer instanton for the splitting between two minima.
//
// A path of P + 1 beads in mass-weighted coordinates q runs from one minimum
// to the other in imaginary time beta hbar, the ends fixed at the minima. The
// instanton is the path that minimises the discretised Euclidean action
//
//   S = sum_j |q_{j+1} - q_j|^2 / (2 dtau) + dtau sum_j' V(q_j),
//
// dtau = beta hbar / P, with half weight on the two ends (trapezoid). The
// splitting follows from the ratio of the off-diagonal to the diagonal
// imaginary-time propagator, both in the same steepest-descent
// approximation:
//
//   delta0 = 2 hbar sqrt(S0 / (2 pi hbar dtau))
//            sqrt(det J_well / det' J) exp(-(S - S_well) / hbar),
//
// with J the Hessian of S over the interior beads, det' leaving out its
// zero mode (the kink's position in imaginary time), S0 the integral of
// |dq/dtau|^2 and J_well the same Hessian with every bead at a minimum.
// Time runs in amu^0.5 Angstrom eV^-0.5, so hbar is kHbar. Where the minimum
// energy path curves, the instanton cuts the corner and follows the
// transverse zero-point energy, which a one-dimensional WKB integral along
// the path cannot.

struct InstantonOptions {
  long beads = 256;             ///< P: segments from one minimum to the other
  double betaHbarOmega = 30.0;  ///< beta hbar omega of the stiffer end along
                                ///< the path; sets the imaginary time
  long maxIterations = 5000;    ///< L-BFGS iterations
  double forceTolerance = 1e-4; ///< largest per-bead |dS/dq| / dtau,
                                ///< eV / (amu^0.5 Angstrom)
  long memory = 20;             ///< L-BFGS correction pairs
};

/// V (eV) and dV/dq (eV / (amu^0.5 Angstrom)) at every point of `q`, all in
/// one call so a potential can spread the beads over its calculators.
using BatchPotential =
    std::function<void(const std::vector<VectorXd> &q, std::vector<double> &v,
                       std::vector<VectorXd> &grad)>;

/// The mass-weighted Hessian d2V/dq2 at interior bead j (1..P-1).
using BeadHessian = std::function<MatrixXd(long j, const VectorXd &q)>;

struct Instanton {
  std::vector<VectorXd> path;   ///< P + 1 beads, ends at the minima
  std::vector<double> energies; ///< V at every bead, eV
  double betaHbar = 0.0;        ///< imaginary time the path spans
  double dtau = 0.0;            ///< betaHbar / P
  double action = 0.0;          ///< (S - S_well) / hbar
  double s0 = 0.0;              ///< integral of |dq/dtau|^2 dtau
  double zeroMode = 0.0;        ///< the eigenvalue det' leaves out
  double delta0 = 0.0;          ///< tunnelling splitting, eV
  double asymmetry = 0.0;       ///< V(end) - V(start), eV
  long iterations = 0;
  bool converged = false;
  /// beta |asymmetry| < 0.1: the propagator ratio reads delta0 only when
  /// the wells lie within a small fraction of kB T of each other.
  bool symmetricEnough = false;
  /// The second smallest eigenvalue of J over the zero mode's: small means
  /// the kink is not isolated in imaginary time and beta hbar is too short.
  double modeSeparation = 0.0;
};

/// omega along the straight line between the minima from the curvature of
/// each well there, the larger of the two; in 1 / time.
double pathOmega(const MatrixXd &hessStart, const MatrixXd &hessEnd,
                 const VectorXd &start, const VectorXd &end);

/// Minimises the action from `guess` (P + 1 beads, ends at the minima, or
/// empty for a tanh kink along the straight line). Evaluates the interior
/// beads once per iteration and line-search step.
Instanton optimizeInstanton(const VectorXd &start, const VectorXd &end,
                            double betaHbar, std::vector<VectorXd> guess,
                            const BatchPotential &potential,
                            const InstantonOptions &options);

/// Fills delta0, zeroMode and modeSeparation from the bead Hessians and the
/// Hessians of the two minima.
void instantonSplitting(Instanton &inst, const BeadHessian &hessian,
                        const MatrixXd &hessStart, const MatrixXd &hessEnd);

} // namespace eonc::tunneling
