//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include <cstdio>
#include <errno.h>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>

#include "IMD.h"

void IMD::cleanMemory(void) {
  system("rm -f imd_eon.out*");
  system("rm -f imd_eon.in.conf");
  return;
}

IMD::~IMD() { cleanMemory(); }

void IMD::force(long N, const double *R, const int *atomicNrs, double *F,
                double *U, double *variance, const double *box) {
  variance = nullptr;
  writeConfIMD(N, R, atomicNrs, box);
  system("imd_eon -p imd_eon.in.param > /dev/null");
  readForceIMD(N, F, U);
  return;
}

void IMD::writeConfIMD(long N, const double *R, const int *atomicNrs,
                       const double *box) {
  int i = 0;
  FILE *confIMD;
  confIMD = fopen("imd_eon.in.conf", "w");

  // boundary box
  fprintf(confIMD, "#X %.8f\t%.8f\t%.8f\n", box[0], box[2], box[3]);
  fprintf(confIMD, "#Y %.8f\t%.8f\t%.8f\n", box[3], box[4], box[5]);
  fprintf(confIMD, "#Z %.8f\t%.8f\t%.8f\n", box[6], box[7], box[8]);

  // header for imd must end with #E
  fprintf(confIMD, "#E\n");

  // positions
  for (i = 0; i < N; i++) {
    fprintf(confIMD, "%i\t%i\t%.5f\t%.19f\t%.19f\t%.19f\n", i, atomicNrs[i] - 1,
            1.0, R[i * 3 + 0], R[i * 3 + 1], R[i * 3 + 2]);
  }
  fclose(confIMD);

  return;
}

void IMD::readForceIMD(long N, double *F, double *U) {
  double junkF;
  char junkChar[256];

  FILE *energyIMD;
  energyIMD = fopen("imd_eon.out.eng", "r");
  double energyPerAtom;
  fgets(junkChar, 256, energyIMD);
  fscanf(energyIMD, "%lf %lf %lf %lf %lf %lf %lf %lf", &junkF, &junkF,
         &energyPerAtom, &junkF, &junkF, &junkF, &junkF, &junkF);
  // energy is given per atom in imd
  *U = N * energyPerAtom;
  fclose(energyIMD);

  FILE *forceIMD;
  forceIMD = fopen("imd_eon.out.00000.wf", "r");
  double forceX;
  double forceY;
  double forceZ;
  double index;
  for (int i = 0; i < N; i++) {
    fscanf(forceIMD, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &index, &junkF,
           &junkF, &junkF, &junkF, &junkF, &forceX, &forceY, &forceZ, &junkF);
    F[int(index) * 3 + 0] = forceX;
    F[int(index) * 3 + 1] = forceY;
    F[int(index) * 3 + 2] = forceZ;
  }
  fclose(forceIMD);
  return;
}
