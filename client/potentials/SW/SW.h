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

#ifndef SW_POTENTIAL
#define SW_POTENTIAL

#include "../../Potential.h"

extern "C" {
/** External function implemented in Fortran to calculate interactions between
atoms using the Stillinger Weber forcefield
@param[in] N           number of atoms
@param[in] R           array to positions of the atoms in Angstrom
@param[out] F           array used to return the forces resulting from
interactions between molecules. Forces are in eV/Angstrom
@param[out] U           pointer to energy in eV
@param[in]  bx, by, bz  pointer to box dimensions in Angstrom
*/
void sw_(const long int *N, const double *R, double *F, double *U,
         const double *bx, const double *by, const double *bz);
}
namespace eonc {
/** SW potential.*/
class SW : public Potential<SW> {
public:
  SW() {}
  void forceImpl(const ForceInput &, ForceOut *) override final;
};
#endif
}
