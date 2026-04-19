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

/** Fortran interface for the Fe-He potential.
@param[in]  N           number of atoms
@param[in]  RX          array of x positions of the atoms in Angstrom
@param[in]  RY          array of y positions of the atoms in Angstrom
@param[in]  RZ          array of z positions of the atoms in Angstrom
@param[in]  ISPEC       array of species identifiers; Fe = 0, He = 1
@param[out] FX          array of x forces between atoms, in eV/Angstrom
@param[out] FY          array of y forces between atoms, in eV/Angstrom
@param[out] FZ          array of z forces between atoms, in eV/Angstrom
@param[out] U           pointer to energy in eV
@param[in]  bx, by, bz  pointer to box dimensions in Angstrom
*/
extern "C" {
void feforce_(const long int *N, const double *RX, const double *RY,
              const double *RZ, const int *ISPEC, double *FX, double *FY,
              double *FZ, double *U, const double *bx, const double *by,
              const double *bz);
}

/// Fe-He interatomic potential (Fortran implementation).
class FeHe final : public Potential {
public:
  explicit FeHe(const Parameters &params)
      : Potential(params) {}
  ~FeHe() override = default;

  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
};
