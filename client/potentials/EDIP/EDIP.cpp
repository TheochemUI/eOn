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

#include "EDIP.h"

void EDIP::cleanMemory(void) { return; }

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void EDIP::force(long N, const double *R, const int *atomicNrs, double *F,
                 double *U, double *variance, const double *box) {
  variance = nullptr;
  edip_(&N, R, F, U, &box[0], &box[4], &box[8]);
  return;
}
