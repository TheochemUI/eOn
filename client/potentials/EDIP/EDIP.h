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

/** External function implemented in Fortran. Calculate interactions between
molecules of water using forcefield EDIP.
@param[in] N           Number of atoms.
@param[in] R           Array to positions of the atoms in Angstrom.
@param[out] F           Array used to return the forces resulting from
interactions between molecules. Forces are in eV/Angstrom.
@param[out] U           Pointer to energy in eV.
@param[in]  bx, by, bz  Pointer to box dimensions in Angstrom.
*/
extern "C" {
void edip_(const long int *N, const double *R, double *F, double *U,
           const double *bx, const double *by, const double *bz);
}

namespace eonc {
class EDIP : public Potential<EDIP> {
public:
  EDIP() {}
  void forceImpl(const ForceInput &, ForceOut *) override final;
};

} // namespace eonc
