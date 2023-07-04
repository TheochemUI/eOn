//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "ExtPot.h"
#include <iostream>
#include <stdio.h>
#include <unistd.h>

void ExtPot::cleanMemory(void) { return; }

ExtPot::~ExtPot() { cleanMemory(); }

void ExtPot::force(long N, const double *R, const int *atomicNrs, double *F,
                   double *U, double *variance, const double *box) {
  variance = nullptr;
  passToSystem(N, R, atomicNrs, box);
  system(eon_extpot_path);
  recieveFromSystem(N, F, U);
  return;
}

void ExtPot::passToSystem(long N, const double *R, const int *atomicNrs,
                          const double *box)
// 'positions' of all particles and box
{
  FILE *out;
  out = fopen("from_eon_to_extpot", "w");

  for (int i = 0; i < 3; i++) {
    fprintf(out, "%.19f\t%.19f\t%.19f\n", box[i * 3 + 0], box[i * 3 + 1],
            box[i * 3 + 2]);
  }

  for (int i = 0; i < N; i++) {
    fprintf(out, "%i\t%.19f\t%.19f\t%.19f\n", atomicNrs[i], R[i * 3 + 0],
            R[i * 3 + 1], R[i * 3 + 2]);
  }
  fclose(out);
  return;
}

void ExtPot::recieveFromSystem(long N, double *F, double *U)
// first line must be the total 'energy', the following lines should be the
// 'forces'
{
  FILE *in;
  in = fopen("from_extpot_to_eon", "r");

  fscanf(in, "%lf", U);

  for (int i = 0; i < N; i++) {
    fscanf(in, "%lf %lf %lf", &F[i * 3 + 0], &F[i * 3 + 1], &F[i * 3 + 2]);
  }

  fclose(in);
  return;
}
