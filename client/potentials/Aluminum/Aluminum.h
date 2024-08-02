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
#include "client/Potential.h"

extern "C" {
/** External function implemented in Fortran
@param[in]	N           number of atoms
@param[in]	R           array to positions of the atoms in Angstrom
@param[out]	F           array used to return the forces between atoms, in
eV/Angstrom
@param[out]	U           pointer to energy in eV
@param[in]  bx, by, bz  pointer to box dimensions in Angstrom
*/
void alforces_(const long int *N, const double *R, double *F, double *U,
            const double *bx, const double *by, const double *bz);
void potinit_();
}

namespace eonc {

/** Aluminum potential.*/
class Aluminum : public Potential<Aluminum> {
public:
  Aluminum() { potinit_(); };
  void forceImpl(const ForceInput &, ForceOut *) override;
};

} // namespace eonc
