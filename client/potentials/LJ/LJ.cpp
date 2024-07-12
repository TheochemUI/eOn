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

#include "LJ.h"
#include "client/potentials/PotHelpers.hpp"
namespace eonc {
// General Functions
void LJ::setParameters(const LJ::Params &p_a) {
  _u0 = p_a.u0;
  _psi = p_a.psi;
  _cutoff_R = p_a.cutoff_R;
  _cutoff_U = calc_cutoffU(p_a);
  return;
}

double LJ::calc_cutoffU(const LJ::Params &p_a) {
  return (4 * p_a.u0 *
          (pow(p_a.psi / p_a.cutoff_R, 12) - pow(p_a.psi / p_a.cutoff_R, 6)));
}

void LJ::forceImpl(const ForceInput &fip, ForceOut *efvd) {
#ifdef EON_CHECKS
  eonc::pot::checkParams(fip);
  eonc::pot::zeroForceOut(fip.nAtoms, efvd);
#endif
  double diffR{0}, diffRX{0}, diffRY{0}, diffRZ{0}, dU{0}, a{0}, b{0};
  for (size_t i = 0; i < fip.nAtoms - 1; i++) {
    for (size_t j = i + 1; j < fip.nAtoms; j++) {
      diffRX = fip.pos[3 * i] - fip.pos[3 * j];
      diffRY = fip.pos[3 * i + 1] - fip.pos[3 * j + 1];
      diffRZ = fip.pos[3 * i + 2] - fip.pos[3 * j + 2];

      // floor = largest integer value less than argument
      diffRX = diffRX - fip.box[0] * floor(diffRX / fip.box[0] + 0.5);
      diffRY = diffRY - fip.box[4] * floor(diffRY / fip.box[4] + 0.5);
      diffRZ = diffRZ - fip.box[8] * floor(diffRZ / fip.box[8] + 0.5);

      diffR = sqrt(diffRX * diffRX + diffRY * diffRY + diffRZ * diffRZ);

      if (diffR < _cutoff_R) {
        // 4u0((psi/r0)^12-(psi/r0)^6)
        a = pow(_psi / diffR, 6);
        b = 4 * _u0 * a;

        efvd->energy += b * (a - 1) - _cutoff_U;

        dU = -6 * b / diffR * (2 * a - 1);
        // F is the negative derivative
        efvd->F[3 * i] = efvd->F[3 * i] - dU * diffRX / diffR;
        efvd->F[3 * i + 1] = efvd->F[3 * i + 1] - dU * diffRY / diffR;
        efvd->F[3 * i + 2] = efvd->F[3 * i + 2] - dU * diffRZ / diffR;

        efvd->F[3 * j] = efvd->F[3 * j] + dU * diffRX / diffR;
        efvd->F[3 * j + 1] = efvd->F[3 * j + 1] + dU * diffRY / diffR;
        efvd->F[3 * j + 2] = efvd->F[3 * j + 2] + dU * diffRZ / diffR;
      }
    }
  }
  return;
}

} // namespace eonc
