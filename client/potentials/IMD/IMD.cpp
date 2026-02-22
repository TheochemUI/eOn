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

#include <cstdio>
#include <errno.h>
#include <filesystem>
#include <iostream>
#include <stdlib.h>

#include "IMD.h"

void IMD::cleanMemory(void) {
  std::error_code ec;
  for (auto &entry : std::filesystem::directory_iterator(".", ec)) {
    auto fname = entry.path().filename().string();
    if (fname.rfind("imd_eon.out", 0) == 0) {
      std::filesystem::remove(entry.path(), ec);
    }
  }
  std::filesystem::remove("imd_eon.in.conf", ec);
  return;
}

IMD::~IMD() { cleanMemory(); }

void IMD::force(long N, const double *R, const int *atomicNrs, double *F,
                double *U, double *variance, const double *box) {
  variance = nullptr;
  writeConfIMD(N, R, atomicNrs, box);
#ifdef _WIN32
  system("imd_eon -p imd_eon.in.param > NUL");
#else
  system("imd_eon -p imd_eon.in.param > /dev/null");
#endif
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
