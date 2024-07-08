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
//-----------------------------------------------------------------------------------
#include "FeHe.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// address to supercell size
void FeHe::force(long N, const double *R, const int *atomicNrs, double *F,
                 double *U, double *variance, const double *box) {
  variance = nullptr;
  assert((int)N > 1);

  long i;
  double *RX = new double[N];
  double *RY = new double[N];
  double *RZ = new double[N];
  double *FX = new double[N];
  double *FY = new double[N];
  double *FZ = new double[N];
  //    long *ISPEC = new long[N];
  int *ISPEC = new int[N];

  for (i = 0; i < N; i++) {
    RX[i] = R[i * 3 + 0];
    RY[i] = R[i * 3 + 1];
    RZ[i] = R[i * 3 + 2];
    assert(atomicNrs[i] == 26 || atomicNrs[i] == 2);
    if (atomicNrs[i] == 26) {
      ISPEC[i] = 0;
    } else {
      ISPEC[i] = 1;
    }
  }

  feforce_(&N, RX, RY, RZ, ISPEC, FX, FY, FZ, U, &box[0], &box[4], &box[8]);

  //    cout <<"forces in c:\n";
  for (i = 0; i < N; i++) {
    F[i * 3] = FX[i];
    F[i * 3 + 1] = FY[i];
    F[i * 3 + 2] = FZ[i];
    //        cout <<F[i*3]<<" "<<F[i*3+1]<<" "<<F[i*3+2]<<endl;
  }
  //    cout <<endl;

  delete[] RX;
  delete[] RY;
  delete[] RZ;
  delete[] FX;
  delete[] FY;
  delete[] FZ;
  delete[] ISPEC;
  return;
}
