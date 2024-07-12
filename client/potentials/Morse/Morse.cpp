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

#include "Morse.h"
#include "../../PotHelpers.hpp"
#include <cassert>
#include <math.h>
namespace eonc {
/** @file
      @brief Morse potential for platinum
      @author Anonymous (possibly A. Pedersen or G. Henkelman ), revision: Jean
   Claude C. Berthet
      @date Unknown, revision: 2010, University of Iceland
      */

/** @class Morse
      @brief Morse potential for platinum.
      Default parameters are for platinum, but parameters for other atoms may be
   specified in setParameters() Morse energy function: @f$
   V(r)=D_e\left[1-e^{-a(r-r_e)}\right]^2 @f$. Equation and notation from
      http://en.wikipedia.org/wiki/Morse_potential .
      */

void Morse::setParameters(const Params &mpar) {
  re_ = mpar.re;
  De_ = mpar.De;
  a_ = mpar.a;
  cutoff_ = mpar.cutoff;
}

void Morse::forceImpl(const ForceInput &params, ForceOut *efvd) {
#ifdef EON_CHECKS
  eonc::pot::checkParams(params);
  eonc::pot::zeroForceOut(params.nAtoms, efvd);
#endif
  double diffR = 0, diffRX{0}, diffRY{0}, diffRZ{0};

  // Initializing en
  for (size_t i = 0; i < params.nAtoms; i++) {
    for (size_t j = i + 1; j < params.nAtoms; j++) {
      diffRX = params.pos[3 * i] - params.pos[3 * j];
      diffRY = params.pos[3 * i + 1] - params.pos[3 * j + 1];
      diffRZ = params.pos[3 * i + 2] - params.pos[3 * j + 2];
      // floor = largest integer value less than argument
      diffRX = diffRX - params.box[0] * floor(diffRX / params.box[0] + 0.5);
      diffRY = diffRY - params.box[4] * floor(diffRY / params.box[4] + 0.5);
      diffRZ = diffRZ - params.box[8] * floor(diffRZ / params.box[8] + 0.5);
      diffR = sqrt(diffRX * diffRX + diffRY * diffRY + diffRZ * diffRZ);
      assert(std::isnormal(diffR));
      if (diffR < cutoff_) {
        double force{0}, energy{0};
        morse(diffR, energy, force);
        efvd->energy += energy;
        efvd->F[3 * i] += force * diffRX / diffR;
        efvd->F[3 * i + 1] += force * diffRY / diffR;
        efvd->F[3 * i + 2] += force * diffRZ / diffR;
        efvd->F[3 * j] -= force * diffRX / diffR;
        efvd->F[3 * j + 1] -= force * diffRY / diffR;
        efvd->F[3 * j + 2] -= force * diffRZ / diffR;
        efvd->energy -= energyCutoff_;
      }
    }
  }
  return;
}

/** Calculate energy and force.
      @param[in] r Distance between two atoms
      @param[out] energy
      @param[out] force -d(energy)/dr
      @see #Morse for equation.
      */
void Morse::morse(double r, double &energy, double &force) {
  double const d = 1 - exp(-a_ * (r - re_));
  energy = De_ * d * d - De_;
  force = 2 * De_ * d * (d - 1) * a_;
}

} // namespace eonc
