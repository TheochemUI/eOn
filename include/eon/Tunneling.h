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

#include <array>
#include <functional>
#include <limits>
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

/// Closed ring of `beads` samples of `path` whose imaginary-time period is
/// `betaHbar`. Bead 0 is the reactant-side turning point and bead N/2 the
/// other; bead N - j repeats bead j. Throws when the path has no barrier,
/// or when the period at the barrier top already exceeds `betaHbar`.
std::vector<VectorXd> ringFromPath(const std::vector<VectorXd> &path,
                                   const std::vector<double> &energies,
                                   double betaHbar, long beads);

/// ln(k), k in 1/time, for the one-dimensional thermal rate along `profile`.
/// `hwReactant` is hbar omega of the reactant well, in eV. Below the barrier
/// the transmission is the WKB factor; above it, the parabolic continuation.
double wkbLogRateAlongPath(const Profile &profile, double beta,
                           double hwReactant);

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

// Ring-polymer instanton for the thermal rate below the crossover
// temperature (Richardson and Althorpe, J. Chem. Phys. 131, 214106 (2009)).
//
// A closed ring of N beads in mass-weighted coordinates q, beta_N = beta / N,
// has the potential
//
//   U_N = sum_j V(q_j) + sum_j |q_{j+1} - q_j|^2 / (2 beta_N^2 hbar^2),
//
// q_N = q_0. Below T_c = hbar omega_b / (2 pi kB), omega_b the imaginary
// frequency at the saddle, the instanton is a first-order saddle of U_N:
// one negative mode, and one zero mode that cycles the beads. The rate is
//
//   k Z_r = (1 / (beta_N hbar)) sqrt(B_N / (2 pi beta_N hbar^2))
//           prod'_k 1 / (beta_N hbar |omega_k|) exp(-beta_N U_N),
//
// B_N = sum_j |q_{j+1} - q_j|^2, omega_k^2 the eigenvalues of the
// mass-weighted Hessian of U_N with the zero mode left out, and Z_r the
// ring-polymer partition function of the harmonic reactant,
//
//   Z_r = exp(-beta V_r) prod_{k=0}^{N-1} prod_m
//         1 / (beta_N hbar sqrt(lambda_m + 4 sin^2(pi k / N) / (beta_N
//         hbar)^2)),
//
// lambda_m the eigenvalues of the reactant's mass-weighted Hessian. Above
// T_c the ring collapses onto the saddle and the formula no longer applies.

/// One unit of time, sqrt(amu Angstrom^2 / eV), in seconds.
inline constexpr double kTimeUnitSeconds = 1.0180505717871193e-14;

/// T_c = hbar omega_b / (2 pi kB) from the mass-weighted Hessian at the
/// saddle, in K; throws when the Hessian has no negative eigenvalue.
double crossoverTemperature(const MatrixXd &hessSaddle);

/// (pi T_c / T) / sin(pi T_c / T). Defined only above T_c. The factor
/// diverges as T approaches T_c from above and tends to 1 at high T.
double parabolicFactor(double temperature, double crossover);

/// ln(k) for classical harmonic transition-state theory, k in 1/time.
/// rigidModes eigenvalues nearest zero are omitted at each Hessian. The
/// saddle's most negative eigenvalue is the barrier mode and leaves the
/// product.
double harmonicTstLogRate(const MatrixXd &hessReactant,
                          const MatrixXd &hessSaddle, double beta,
                          double barrier, long rigidModes);

/// ln(k) for quantum harmonic transition-state theory, k in 1/time:
/// (1 / (2 pi beta hbar)) prod_r 2 sinh(beta hbar omega_r / 2) /
/// prod'_s 2 sinh(beta hbar omega_s / 2) exp(-beta barrier), the
/// zero-point and quantised partition functions of every bound mode; the
/// saddle's unstable mode leaves the product. Times the parabolic factor it
/// is the N -> infinity ring-polymer rate above T_c, so it joins the
/// instanton rate at the crossover; its high-temperature limit is
/// harmonicTstLogRate.
double quantumHarmonicTstLogRate(const MatrixXd &hessReactant,
                                 const MatrixXd &hessSaddle, double beta,
                                 double barrier, long rigidModes);

struct RateInstantonOptions {
  long beads = 32;              ///< N, beads on the ring
  long maxIterations = 1000;    ///< translation steps
  double forceTolerance = 1e-3; ///< largest per-bead |dU_N/dq|,
                                ///< eV / (amu^0.5 Angstrom)
  long lanczosFirst = 30;       ///< Lanczos steps for the first minimum mode
  long lanczosRestart = 6;      ///< Lanczos steps from the previous mode
  double lanczosStep = 1e-4;    ///< finite-difference step, amu^0.5 Angstrom
  double maxStep = 0.05;        ///< largest bead move per step,
                                ///< amu^0.5 Angstrom
  long memory = 10;             ///< L-BFGS correction pairs
  /// Even N, when the guess already matches under j -> N - j: evaluate
  /// the potential from one turning point to the other and copy it.
  /// An odd N, or a guess without that symmetry, evaluates every bead.
  bool halfRing = true;
  /// A half ring that has converged, or that is stationary at the wrong
  /// index, is probed for an unstable mode odd under j -> N - j (two
  /// copies of the instanton on one ring), at about lanczosFirst * N
  /// gradient calls; the search then finishes on the whole ring. A cooling
  /// schedule probes only its last temperature.
  bool checkOddSector = true;
  double energyShift = 0.0; ///< subtracted from every bead potential, eV
  /// Active coordinates at or below this take the Newton step. Zero keeps
  /// minimum-mode following. The step solves through the block chain, so
  /// the default admits every size a batch potential can evaluate.
  long newtonLimit = std::numeric_limits<long>::max();
  /// Where the bead Hessian blocks start: "saddle" copies the saddle's
  /// Hessian to every bead and lets the Bofill update carry it, at no force
  /// calls; "finite_difference" takes 2 f gradient calls per bead first.
  std::string initialHessians = "saddle";
  /// Rigid motions of the ring: a translation, or one rotation about the
  /// ring's centre of mass, applied to every bead alike leaves U_N
  /// unchanged. The search rebuilds those directions from the current beads
  /// at every step (the rigid quotient) and keeps them out of the step and
  /// out of the classification, as it does the imaginary-time cycle. Empty
  /// masses switch it off (atoms fixed). rigidSqrtMasses holds sqrt(m) per
  /// atom, rigidReference the Cartesian positions q is measured from (3 per
  /// atom), rigidRotations which rotations are free (a free cluster).
  std::vector<double> rigidSqrtMasses;
  VectorXd rigidReference;
  std::array<bool, 3> rigidRotations{{false, false, false}};
};

/// Spectrum of a closed ring's Hessian without forming it.
struct RingSpectrum {
  /// ln |det' J|: the product over every eigenvalue but the one along tau.
  double logDetPrime = 0.0;
  /// Eigenvalues below zero, the one along tau left out.
  long negativeModes = 0;
  /// tau . J tau, the eigenvalue the prime leaves out; small when the ring
  /// is a converged instanton.
  double zeroEigenvalue = 0.0;
};

/// The ring Hessian of bead Hessians `beadHessians` (d2V/dq2 at each of the
/// N beads) and spring constant c, with the normalised direction `tau`
/// (N beads) projected out through the determinant lemma
/// det(J + tau tau^T) = det' J when J tau = 0. Block LU of the open chain
/// plus a low-rank correction for the closure and tau, O(N f^3).
RingSpectrum ringSpectrum(const std::vector<MatrixXd> &beadHessians, double c,
                          const std::vector<VectorXd> &tau);

struct RateInstanton {
  std::vector<VectorXd> beads;     ///< N beads, q_N = q_0 implied
  std::vector<double> energies;    ///< V at every bead, eV
  double beta = 0.0;               ///< 1 / (kB T), 1 / eV
  double betaN = 0.0;              ///< beta / N
  double temperature = 0.0;        ///< K
  double crossover = 0.0;          ///< T_c, K
  double ringPotential = 0.0;      ///< U_N, eV
  double bN = 0.0;                 ///< sum_j |q_{j+1} - q_j|^2, amu Angstrom^2
  double negativeEigenvalue = 0.0; ///< of the ring Hessian, 1 / time^2
  double zeroEigenvalue = 0.0;     ///< the eigenvalue left out
  long negativeModes = 0;          ///< eigenvalues below the zero mode
  long iterations = 0;
  bool converged = false;
  double logRateTimesZr = 0.0; ///< ln(k Z_r), k in 1 / time
  double logZr = 0.0;          ///< ln Z_r
  double logRate = 0.0;        ///< ln k, k in 1 / time
  double rate = 0.0;           ///< k in 1 / s
  /// -kB T ln(2 pi hbar beta k): the barrier an Eyring rate would need, eV.
  double effectiveBarrier = 0.0;
  /// Classical harmonic transition-state theory at the same T, 1 / s,
  /// when the saddle Hessian was given; its logarithm (k in 1 / time) does
  /// not underflow.
  double classicalRate = 0.0;
  double classicalLogRate = 0.0;
};

/// Finds the rate instanton at inverse temperature `beta` (1 / eV).
/// `guess` holds N beads, or is empty for a ring stretched along the saddle's
/// unstable mode to where V has dropped by (1 - T / T_c) of the lower of the
/// two barriers. At or below `newtonLimit` active coordinates the step is an
/// index-1 Newton step on a Bofill Hessian, and below 0.75 T_c an empty
/// guess cools from 0.85 T_c. Beyond that limit the step is minimum-mode
/// following. An even ring whose beads match under j -> N - j is evaluated
/// from one turning point to the other and mirrored. `saddle` and
/// `hessSaddle` are mass-weighted.
RateInstanton optimizeRateInstanton(const VectorXd &saddle,
                                    const MatrixXd &hessSaddle, double beta,
                                    std::vector<VectorXd> guess,
                                    const BatchPotential &potential,
                                    const RateInstantonOptions &options);

/// Where a ring sits against the saddle it was seeded from, in
/// mass-weighted coordinates. s = (q - saddle) . mode along the unstable
/// mode; the turning points are the beads of least and greatest s.
struct RingChannel {
  double sMin = 0.0; ///< least s over the beads, amu^0.5 Angstrom
  double sMax = 0.0; ///< greatest s over the beads
  /// |cos| of the angle between the chord joining the turning points and
  /// the unstable mode.
  double chordOverlap = 0.0;
  /// Largest distance from the saddle, across the mode, of a point where
  /// the ring crosses the dividing plane s = 0.
  double crossingOffset = 0.0;
  /// Beads on both sides of the dividing plane, the chord within 60 degrees
  /// of the mode, and every crossing within sMax - sMin of the saddle. A
  /// ring around a neighbouring saddle can still straddle the plane, but it
  /// crosses it far from this saddle.
  bool belongs = false;
};

/// The channel test above for a ring of beads, a saddle and its unstable
/// mode (normalised inside).
RingChannel ringChannel(const std::vector<VectorXd> &beads,
                        const VectorXd &saddle, const VectorXd &mode);

/// The mass-weighted Hessian d2V/dq2 at ring bead j (0..N-1).
using RingBeadHessian = std::function<MatrixXd(long j, const VectorXd &q)>;

/// The atoms behind a ring's mass-weighted coordinates, for its rigid
/// motions: sqrt(m) per atom, the Cartesian positions q is measured from
/// (3 per atom), and which rotations are free. Empty masses: the rigid
/// directions come from the reactant Hessian instead.
struct RingRigidBodies {
  std::vector<double> sqrtMasses;
  VectorXd reference;
  std::array<bool, 3> rotations{{false, false, false}};
};

/// Fills the rate from the bead Hessians, the reactant minimum's Hessian and
/// energy, and optionally the saddle's Hessian and energy for the classical
/// comparison (pass an empty matrix to skip it). rigidModes is the count of
/// rigid-body zero modes to omit: the translations, plus a rotation only when
/// the reactant Hessian leaves it null (a free cluster has them, a crystal
/// does not, and an atom held fixed has none). They leave the centroid
/// factors, so the rotational and translational partition functions of
/// reactant and instanton cancel; for rotations that is an approximation,
/// since the ring's moments of inertia are not the reactant's. On the ring
/// the omitted directions are its null vectors: with `rigidBodies` given,
/// the translations and rotations of the beads themselves about the ring's
/// centre of mass (a rotation moves each bead differently), otherwise the
/// reactant Hessian's null vectors copied to every bead, which are exact
/// only for a coordinate every bead shares. The bead Hessians must then
/// keep their rotational part: projecting a rotation out of a bead that is
/// not stationary removes the curvature that balances its springs.
///
/// Up to denseLimit ring degrees of freedom, and for any negative limit,
/// the determinant, the negative-mode count and the lowest eigenvalue come
/// from the dense eigenvalues of the lifted ring Hessian, and a block-chain
/// determinant that disagrees with them throws. Beyond that, and for a
/// limit of 0, the block chain gives the determinant and the inertia, and a
/// Lanczos run on ring products the lowest eigenvalue.
void instantonRate(RateInstanton &inst, const RingBeadHessian &hessian,
                   const MatrixXd &hessReactant, double vReactant,
                   const MatrixXd &hessSaddle = MatrixXd(),
                   double vSaddle = 0.0, long rigidModes = 0,
                   long denseLimit = 4096,
                   const RingRigidBodies &rigidBodies = {});

/// log|det| of the cyclic block-tridiagonal ring Hessian. Each diag[j]
/// already contains the bead Hessian plus 2 c I, and the neighbour coupling
/// is -c I, including the corner that closes the ring. A singular ring
/// returns -infinity. Throws when the open chain is singular.
double cyclicRingLogAbsDet(double c, const std::vector<MatrixXd> &diag);

/// Solves that same cyclic ring Hessian. Throws when the ring is singular
/// or the right-hand side does not match the blocks.
std::vector<VectorXd> cyclicRingSolve(double c,
                                      const std::vector<MatrixXd> &diag,
                                      const std::vector<VectorXd> &rhs);

} // namespace eonc::tunneling
