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

/** @file
      @brief Morse potential for platinum
      @author Anonymous (possibly A. Pedersen or G. Henkelman), revision: Jean
   Claude C. Berthet
      @date Unknown, revision: 2010, University of Iceland
      */
#include "../../Potential.h"
#include <cmath>

/// Morse potential, default parameters are for Pt.
/// V(r) = De * [1 - exp(-a*(r - re))]^2
class Morse final : public Potential {
public:
  explicit Morse(const Parameters &params)
      : Potential(params), De_{0.7102}, a_{1.6047}, re_{2.8970}, cutoff_{9.5} {
    setParameters(De_, a_, re_, cutoff_);
  }

  /// Parameters: De in eV, a in 1/Angstrom, re in Angstrom, cutoff in Angstrom.
  void setParameters(double De, double a, double re, double cutoff);

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  /// Compute Morse energy and force magnitude for distance r.
  /// @param[in]  r      distance between two atoms
  /// @param[out] energy Morse energy at r
  /// @param[out] force  -d(energy)/dr (force magnitude)
  void morse(double r, double &energy, double &force) const;

  double De_{0.0};
  double a_{0.0};
  double re_{0.0};
  double cutoff_{0.0};
  double energyCutoff_{0.0};
};
