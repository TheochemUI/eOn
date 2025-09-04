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
atoms using the Stillinger Weber forcefield
@param[in]	N           number of atoms
@param[in]	R           array to positions of the atoms in Angstrom
@param[out]	F           array used to return the forces resulting from
interactions between molecules. Forces are in eV/Angstrom
@param[out]	U           pointer to energy in eV
@param[in]  bx, by, bz  pointer to box dimensions in Angstrom
*/
extern "C" {
void sw_(const long int *N, const double *R, double *F, double *U,
         const double *bx, const double *by, const double *bz);
}

/** SW potential.*/
class SW : public Potential {
private:
  std::shared_ptr<Parameters> parameters;

public:
  // Functions
  // constructor
  SW(std::shared_ptr<Parameters> p)
      : Potential(p),
        parameters{p} {}

  // To satisfy interface
  void initialize(void);
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);
};
