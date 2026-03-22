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

/** Fortran interface for the Stillinger-Weber potential.
@param[in]  N           number of atoms
@param[in]  R           array of positions in Angstrom
@param[out] F           array of forces in eV/Angstrom
@param[out] U           pointer to energy in eV
@param[in]  bx, by, bz pointer to box dimensions in Angstrom
*/
extern "C" {
void sw_(const long int *N, const double *R, double *F, double *U,
         const double *bx, const double *by, const double *bz);
}

/// Stillinger-Weber potential for Si (Fortran implementation).
class SW final : public Potential {
public:
  explicit SW(const Parameters &p) : Potential(p) {}
  ~SW() override = default;

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
};
