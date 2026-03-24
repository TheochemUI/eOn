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

#include "../../Potential.h"

/** Fortran interface for the Tersoff potential.
@param[in]  N           Number of atoms
@param[in]  R           Array of positions in Angstrom
@param[out] F           Array of forces in eV/Angstrom
@param[out] U           Pointer to energy in eV
@param[in]  bx, by, bz Pointer to box dimensions in Angstrom
*/
extern "C" {
void tersoff_(const long int *N, const double *R, double *F, double *U,
              const double *bx, const double *by, const double *bz);
}

/// Tersoff potential for Si (Fortran implementation).
class Tersoff final : public Potential {
public:
  explicit Tersoff(const Parameters &p)
      : Potential(p) {}
  ~Tersoff() override = default;

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
};
