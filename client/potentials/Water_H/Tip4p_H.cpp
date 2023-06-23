//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

// the add H must be the last atom
#include "Tip4p_H.h"

void Tip4p_H::cleanMemory(void) { return; }

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// address to supercell size
void Tip4p_H::force(long N, const double *R, const int *atomicNrs, double *F,
                    double *U, double *variance, const double *box) {
  variance = nullptr;

  double diffR[3], force[3];
  double energy;
  double diffRX, diffRY, diffRZ;
  double diffR_O = 0, diffRX_O, diffRY_O, diffRZ_O;
  double diffR_H1 = 0, diffRX_H1, diffRY_H1, diffRZ_H1;
  double diffR_H2 = 0, diffRX_H2, diffRY_H2, diffRZ_H2;

  int nMolecules = (N - 1) / 3;
  int indexStartO = 6 * nMolecules;
  int indexAddH = 3 * (N - 1);

  *U = 0;
  energy = 0;

  for (int i = 0; i < N; i++) {
    F[3 * i] = 0;
    F[3 * i + 1] = 0;
    F[3 * i + 2] = 0;
  }
  double posX_H, posY_H, posZ_H;
  posX_H = R[indexAddH];
  posY_H = R[indexAddH + 1];
  posZ_H = R[indexAddH + 2];
  // initializing end

  // remember that there is an additional H atom
  // water - water interaction
  tip4p_pot->force(N - 1, R, atomicNrs, F, U, box);

  for (int i = 0; i < nMolecules; i++) {
    // hydrogen atoms are located together in pairs of 2 (6 coordinates)
    // h1 atom
    diffRX = R[6 * i] - posX_H;
    diffRY = R[6 * i + 1] - posY_H;
    diffRZ = R[6 * i + 2] - posZ_H;

    // floor = largest integer value less than argument

    diffRX_H1 = diffRX - box[0] * floor(diffRX / box[0] + 0.5);
    diffRY_H1 = diffRY - box[4] * floor(diffRY / box[4] + 0.5);
    diffRZ_H1 = diffRZ - box[8] * floor(diffRZ / box[8] + 0.5);
    diffR_H1 = sqrt(diffRX_H1 * diffRX_H1 + diffRY_H1 * diffRY_H1 +
                    diffRZ_H1 * diffRZ_H1);

    // h2 atom
    diffRX = R[6 * i + 3] - posX_H;
    diffRY = R[6 * i + 4] - posY_H;
    diffRZ = R[6 * i + 5] - posZ_H;
    diffRX_H2 = diffRX - box[0] * floor(diffRX / box[0] + 0.5);
    diffRY_H2 = diffRY - box[4] * floor(diffRY / box[4] + 0.5);
    diffRZ_H2 = diffRZ - box[8] * floor(diffRZ / box[8] + 0.5);
    diffR_H2 = sqrt(diffRX_H2 * diffRX_H2 + diffRY_H2 * diffRY_H2 +
                    diffRZ_H2 * diffRZ_H2);

    // oxygen atoms are in the end
    diffRX = R[3 * i + indexStartO] - posX_H;
    diffRY = R[3 * i + 1 + indexStartO] - posY_H;
    diffRZ = R[3 * i + 2 + indexStartO] - posZ_H;
    diffRX_O = diffRX - box[0] * floor(diffRX / box[0] + 0.5);
    diffRY_O = diffRY - box[4] * floor(diffRY / box[4] + 0.5);
    diffRZ_O = diffRZ - box[8] * floor(diffRZ / box[8] + 0.5);
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

    *U += energy / 96.485336; // kJ/mol to eV
    force[0] /= 96.485336;
    force[1] /= 96.485336;
    force[2] /= 96.485336;

    // oxygen atom
    F[3 * i + indexStartO] -= force[0] * diffRX_O / diffR_O;
    F[3 * i + 1 + indexStartO] -= force[0] * diffRY_O / diffR_O;
    F[3 * i + 2 + indexStartO] -= force[0] * diffRZ_O / diffR_O;
    F[indexAddH] += force[0] * diffRX_O / diffR_O;
    F[indexAddH + 1] += force[0] * diffRY_O / diffR_O;
    F[indexAddH + 2] += force[0] * diffRZ_O / diffR_O;

    // h1 atom
    F[6 * i] -= force[1] * diffRX_H1 / diffR_H1;
    F[6 * i + 1] -= force[1] * diffRY_H1 / diffR_H1;
    F[6 * i + 2] -= force[1] * diffRZ_H1 / diffR_H1;
    F[indexAddH] += force[1] * diffRX_H1 / diffR_H1;
    F[indexAddH + 1] += force[1] * diffRY_H1 / diffR_H1;
    F[indexAddH + 2] += force[1] * diffRZ_H1 / diffR_H1;

    // h2 atom
    F[6 * i + 3] -= force[2] * diffRX_H2 / diffR_H2;
    F[6 * i + 4] -= force[2] * diffRY_H2 / diffR_H2;
    F[6 * i + 5] -= force[2] * diffRZ_H2 / diffR_H2;
    F[indexAddH] += force[2] * diffRX_H2 / diffR_H2;
    F[indexAddH + 1] += force[2] * diffRY_H2 / diffR_H2;
    F[indexAddH + 2] += force[2] * diffRZ_H2 / diffR_H2;
  }

  return;
}
