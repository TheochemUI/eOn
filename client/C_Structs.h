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
  /**
  @param[in] nAtoms           Number of atoms.
  @param[in] pos              Pointer to an array of positions Angstrom.
  @param[in] atmnrs           Pointer to an array of atomic numbers.
  @param[in] box              Pointer to supercell dimensions in Angstrom.
  */
  const size_t nAtoms;
  const double *pos;
  const size_t *atmnrs;
  const double *box;
} ForceInput;

typedef struct {
  /**
  @param[out] F           Array used to return the forces resulting from
  interactions between molecules. Forces are in eV/Angstrom.
  @param[out] energy      Energy in eV.
  @param[out] variance    Variance, 0 when not needed.
   */
  double *F;
  double energy;
  double variance;
} ForceOut;

#ifdef __cplusplus
}
#endif
