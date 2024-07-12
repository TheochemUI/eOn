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

#pragma once

#include <stddef.h>

// These live outside the namespace as well
// Remember, they are C-style!!

#ifdef __cplusplus
extern "C" {
#endif

// These are the only structs without member initialization since they are
// C-style
typedef struct {
  // pointer to number of atoms, pointer to array of positions
  // address to supercell size
  const size_t nAtoms;
  const double *pos;
  const size_t *atmnrs;
  const double *box;
} ForceInput;

typedef struct {
  // pointer to array of forces
  double *F;
  // Internal energy
  double energy;
  // Variance here is 0 when not needed and that's OK
  double variance;
} ForceOut;

#ifdef __cplusplus
}
#endif
