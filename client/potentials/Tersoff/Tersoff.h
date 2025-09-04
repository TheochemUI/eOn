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

/** External function implemented in Fortran to calculate interactions between
atoms using the Tersoff forcefield
@param[in]	N           Number of atoms
@param[in]	R           Array to positions of the atoms in Angstrom
@param[out]	F           Array used to return the forces between atoms, in
eV/Angstrom
@param[out]	U           Pointer to energy in eV
@param[in]  bx, by, bz  Pointer to box dimensions in Angstrom
*/
extern "C" {
void tersoff_(const long int *N, const double *R, double *F, double *U,
              const double *bx, const double *by, const double *bz);
}

/** Tersoff potential */
class Tersoff : public Potential {

public:
  // Functions
  // constructor
  Tersoff(std::shared_ptr<Parameters> p)
      : Potential(p) {}

  // To satisfy interface
  void initialize(void);
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);
};
