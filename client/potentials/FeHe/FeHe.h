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

/** External function implemented in Fortran
@param[in]	N           number of atoms
@param[in]	RX          array of x positions of the atoms in Angstrom
@param[in]	RY          array of x positions of the atoms in Angstrom
@param[in]	RZ          array of x positions of the atoms in Angstrom
@param[in]	ISPEC       array of species itentifiers; Fe = 0, He = 1
@param[out]	FX          array used to return the x forces between atoms, in
eV/Angstrom
@param[out]	FY          array used to return the y forces between atoms, in
eV/Angstrom
@param[out]	FZ          array used to return the z forces between atoms, in
eV/Angstrom
@param[out]	U           pointer to energy in eV
@param[in]  bx, by, bz  pointer to box dimensions in Angstrom
*/

extern "C" {
void feforce_(const long int *N, const double *RX, const double *RY,
              const double *RZ, const int *ISPEC, double *FX, double *FY,
              double *FZ, double *U, const double *bx, const double *by,
              const double *bz);
}
/* No potential initialization
extern "C"
{
    void potinit_();
}
*/
/** FeHe potential.*/
class FeHe : public Potential {
public:
  FeHe(std::shared_ptr<Parameters> params)
      : Potential(params) {}
  ~FeHe(void) {};
  // To satisfy interface
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
};
