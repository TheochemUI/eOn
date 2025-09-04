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
