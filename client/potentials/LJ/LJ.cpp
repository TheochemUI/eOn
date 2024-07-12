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
#include "../../PotHelpers.hpp"
namespace eonc {
// General Functions
void LJ::setParameters(const Params &ljp) {
  u0 = ljp.u0;
  psi = ljp.psi;
  cutoff_R = ljp.cutoff_R;
  cutoff_U = calc_cutoffU(ljp);
  return;
}

double LJ::calc_cutoffU(const LJ::Params &p) {
  return (4 * u0 * (pow(psi / cutoff_R, 12) - pow(psi / cutoff_R, 6)));
}

void LJ::forceImpl(const ForceInput &params, ForceOut *efvd) {
#ifdef EON_CHECKS
  eonc::pot::zeroForceOut(params.nAtoms, efvd);
#endif
  double diffR{0}, diffRX{0}, diffRY{0}, diffRZ{0}, dU{0}, a{0}, b{0};
  for (size_t i = 0; i < params.nAtoms - 1; i++) {
    for (size_t j = i + 1; j < params.nAtoms; j++) {
      diffRX = params.pos[3 * i] - params.pos[3 * j];
      diffRY = params.pos[3 * i + 1] - params.pos[3 * j + 1];
      diffRZ = params.pos[3 * i + 2] - params.pos[3 * j + 2];

      // floor = largest integer value less than argument
      diffRX = diffRX - params.box[0] * floor(diffRX / params.box[0] + 0.5);
      diffRY = diffRY - params.box[4] * floor(diffRY / params.box[4] + 0.5);
      diffRZ = diffRZ - params.box[8] * floor(diffRZ / params.box[8] + 0.5);

      diffR = sqrt(diffRX * diffRX + diffRY * diffRY + diffRZ * diffRZ);

      if (diffR < cutoff_R) {
        // 4u0((psi/r0)^12-(psi/r0)^6)
        a = pow(psi / diffR, 6);
        b = 4 * u0 * a;

        efvd->energy = efvd->energy + b * (a - 1) - cutoff_U;

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
