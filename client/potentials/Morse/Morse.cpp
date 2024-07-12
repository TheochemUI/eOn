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
#include "client/potentials/PotHelpers.hpp"
#include <cmath>
namespace eonc {
/** @class Morse
      @brief Morse potential for platinum.
      Default parameters are for platinum, but parameters for other atoms may be
   specified in setParameters() Morse energy function: @f$
   V(r)=D_e\left[1-e^{-a(r-r_e)}\right]^2 @f$. Equation and notation from
      http://en.wikipedia.org/wiki/Morse_potential .
      */

void Morse::setParameters(const Morse::Params &p_a) {
  _re = p_a.re;
  _De = p_a.De;
  _a = p_a.a;
  _cutoff = p_a.cutoff;
}

void Morse::forceImpl(const ForceInput &fip, ForceOut *efvd) {
#ifdef EON_CHECKS
  eonc::pot::checkParams(fip);
  eonc::pot::zeroForceOut(fip.nAtoms, efvd);
#endif
  double diffR = 0, diffRX{0}, diffRY{0}, diffRZ{0};

  // Initializing en
  for (size_t i = 0; i < fip.nAtoms; i++) {
    for (size_t j = i + 1; j < fip.nAtoms; j++) {
      diffRX = fip.pos[3 * i] - fip.pos[3 * j];
      diffRY = fip.pos[3 * i + 1] - fip.pos[3 * j + 1];
      diffRZ = fip.pos[3 * i + 2] - fip.pos[3 * j + 2];
      // floor = largest integer value less than argument
      diffRX = diffRX - fip.box[0] * floor(diffRX / fip.box[0] + 0.5);
      diffRY = diffRY - fip.box[4] * floor(diffRY / fip.box[4] + 0.5);
      diffRZ = diffRZ - fip.box[8] * floor(diffRZ / fip.box[8] + 0.5);
      diffR = sqrt(diffRX * diffRX + diffRY * diffRY + diffRZ * diffRZ);
      assert(std::isnormal(diffR));
      if (diffR < _cutoff) {
        double force{0}, energy{0};
        morse(diffR, energy, force);
        efvd->energy += energy;
        efvd->F[3 * i] += force * diffRX / diffR;
        efvd->F[3 * i + 1] += force * diffRY / diffR;
        efvd->F[3 * i + 2] += force * diffRZ / diffR;
        efvd->F[3 * j] -= force * diffRX / diffR;
        efvd->F[3 * j + 1] -= force * diffRY / diffR;
        efvd->F[3 * j + 2] -= force * diffRZ / diffR;
        efvd->energy -= _energyCutoff;
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
  double const d = 1 - exp(-_a * (r - _re));
  energy = _De * d * d - _De;
  force = 2 * _De * d * (d - 1) * _a;
}

} // namespace eonc
