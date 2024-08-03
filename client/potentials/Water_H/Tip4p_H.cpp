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

// the add H must be the last atom
#include "Tip4p_H.h"
namespace eonc {
void Tip4p_H::forceImpl(const ForceInput &fip, ForceOut *efvd) {
#ifdef EON_CHECKS
  eonc::pot::checkParams(fip);
  eonc::pot::zeroForceOut(fip.nAtoms, efvd);
#endif
  const long int N = fip.nAtoms;

  double diffR[3], force[3];
  double energy{0};
  double diffRX{0}, diffRY{0}, diffRZ{0};
  double diffR_O{0}, diffRX_O{0}, diffRY_O{0}, diffRZ_O{0};
  double diffR_H1{0}, diffRX_H1{0}, diffRY_H1{0}, diffRZ_H1{0};
  double diffR_H2{0}, diffRX_H2{0}, diffRY_H2{0}, diffRZ_H2{0};

  int nMolecules = (N - 1) / 3;
  int indexStartO = 6 * nMolecules;
  int indexAddH = 3 * (N - 1);

  double posX_H{0}, posY_H{0}, posZ_H{0};
  posX_H = fip.pos[indexAddH];
  posY_H = fip.pos[indexAddH + 1];
  posZ_H = fip.pos[indexAddH + 2];
  // initializing end

  // remember that there is an additional H atom
  // water - water interaction
  ForceInput fwatwat{static_cast<size_t>(N - 1), fip.pos, fip.atmnrs, fip.box};
  tip4p_pot->forceImpl(fwatwat, efvd);

  for (int i = 0; i < nMolecules; i++) {
    // hydrogen atoms are located together in pairs of 2 (6 coordinates)
    // h1 atom
    diffRX = fip.pos[6 * i] - posX_H;
    diffRY = fip.pos[6 * i + 1] - posY_H;
    diffRZ = fip.pos[6 * i + 2] - posZ_H;

    // floor = largest integer value less than argument

    diffRX_H1 = diffRX - fip.box[0] * floor(diffRX / fip.box[0] + 0.5);
    diffRY_H1 = diffRY - fip.box[4] * floor(diffRY / fip.box[4] + 0.5);
    diffRZ_H1 = diffRZ - fip.box[8] * floor(diffRZ / fip.box[8] + 0.5);
    diffR_H1 = sqrt(diffRX_H1 * diffRX_H1 + diffRY_H1 * diffRY_H1 +
                    diffRZ_H1 * diffRZ_H1);

    // h2 atom
    diffRX = fip.pos[6 * i + 3] - posX_H;
    diffRY = fip.pos[6 * i + 4] - posY_H;
    diffRZ = fip.pos[6 * i + 5] - posZ_H;
    diffRX_H2 = diffRX - fip.box[0] * floor(diffRX / fip.box[0] + 0.5);
    diffRY_H2 = diffRY - fip.box[4] * floor(diffRY / fip.box[4] + 0.5);
    diffRZ_H2 = diffRZ - fip.box[8] * floor(diffRZ / fip.box[8] + 0.5);
    diffR_H2 = sqrt(diffRX_H2 * diffRX_H2 + diffRY_H2 * diffRY_H2 +
                    diffRZ_H2 * diffRZ_H2);

    // oxygen atoms are in the end
    diffRX = fip.pos[3 * i + indexStartO] - posX_H;
    diffRY = fip.pos[3 * i + 1 + indexStartO] - posY_H;
    diffRZ = fip.pos[3 * i + 2 + indexStartO] - posZ_H;
    diffRX_O = diffRX - fip.box[0] * floor(diffRX / fip.box[0] + 0.5);
    diffRY_O = diffRY - fip.box[4] * floor(diffRY / fip.box[4] + 0.5);
    diffRZ_O = diffRZ - fip.box[8] * floor(diffRZ / fip.box[8] + 0.5);
    diffR_O =
        sqrt(diffRX_O * diffRX_O + diffRY_O * diffRY_O + diffRZ_O * diffRZ_O);

    diffR[0] = diffR_O;
    diffR[1] = diffR_H1;
    diffR[2] = diffR_H2;
    force[0] = 0;
    force[1] = 0;
    force[2] = 0;
    energy = 0;

    // h - water interaction
    poth2oh_(diffR, &energy, force);

    efvd->energy += energy / 96.485336; // kJ/mol to eV
    force[0] /= 96.485336;
    force[1] /= 96.485336;
    force[2] /= 96.485336;

    // oxygen atom
    efvd->F[3 * i + indexStartO] -= force[0] * diffRX_O / diffR_O;
    efvd->F[3 * i + 1 + indexStartO] -= force[0] * diffRY_O / diffR_O;
    efvd->F[3 * i + 2 + indexStartO] -= force[0] * diffRZ_O / diffR_O;
    efvd->F[indexAddH] += force[0] * diffRX_O / diffR_O;
    efvd->F[indexAddH + 1] += force[0] * diffRY_O / diffR_O;
    efvd->F[indexAddH + 2] += force[0] * diffRZ_O / diffR_O;

    // h1 atom
    efvd->F[6 * i] -= force[1] * diffRX_H1 / diffR_H1;
    efvd->F[6 * i + 1] -= force[1] * diffRY_H1 / diffR_H1;
    efvd->F[6 * i + 2] -= force[1] * diffRZ_H1 / diffR_H1;
    efvd->F[indexAddH] += force[1] * diffRX_H1 / diffR_H1;
    efvd->F[indexAddH + 1] += force[1] * diffRY_H1 / diffR_H1;
    efvd->F[indexAddH + 2] += force[1] * diffRZ_H1 / diffR_H1;

    // h2 atom
    efvd->F[6 * i + 3] -= force[2] * diffRX_H2 / diffR_H2;
    efvd->F[6 * i + 4] -= force[2] * diffRY_H2 / diffR_H2;
    efvd->F[6 * i + 5] -= force[2] * diffRZ_H2 / diffR_H2;
    efvd->F[indexAddH] += force[2] * diffRX_H2 / diffR_H2;
    efvd->F[indexAddH + 1] += force[2] * diffRY_H2 / diffR_H2;
    efvd->F[indexAddH + 2] += force[2] * diffRZ_H2 / diffR_H2;
  }

  return;
}

} // namespace eonc
