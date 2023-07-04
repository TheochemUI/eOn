//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "NewPot.h"

void NewPot::cleanMemory(void) { return; }

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void NewPot::force(long N, const double *R, const int *atomicNrs, double *F,
                   double *U, double *variance, const double *box) {
  variance = nullptr;

  for (int i = 0; i < N; i++) {
    F[3 * i] = fake1;
    F[3 * i + 1] = fake1;
    F[3 * i + 2] = fake1;
  }

  *U = fake2;
  return;
}
