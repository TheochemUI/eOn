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

#include <memory>
#include <string>
#include <vector>

namespace eonc::tunneling {

/// hbar in eV^0.5 amu^0.5 Angstrom: 1.054571817e-34 J s over
/// sqrt(1.602176634e-19 J * 1.66053906660e-27 kg) * 1e-10 m.
inline constexpr double kHbar = 0.06465415130134121;

/// Boltzmann constant in eV / K.
inline constexpr double kBoltzmann = 8.617333262e-5;

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

} // namespace eonc::tunneling
