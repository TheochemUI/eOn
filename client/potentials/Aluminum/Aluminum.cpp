//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//
//-----------------------------------------------------------------------------------
#include "Aluminum.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// address to supercell size
void Aluminum::force(long N, const double *R, const int *atomicNrs, double *F,
                     double *U, double *variance, const double *box) {
  assert((int)N > 1);
  force_(&N, R, F, U, &box[0], &box[4], &box[8]);
  variance = nullptr;
  return;
}
