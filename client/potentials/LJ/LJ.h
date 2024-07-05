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
// #include "../../system_unit.h" // unit converters
#include "../../Potential.h"

namespace eonc::def {
struct LJParams {
  double u0;
  double cutoff;
  double psi;
  LJParams()
      : u0{1.0},
        cutoff{15.0},
        psi{1.0} {}
};
} // namespace eonc::def

/** Lennard Jones potential.*/
class LJ : public Potential {
private:
  double u0;
  double cuttOffR;
  double psi;
  double cuttOffU;

public:
  LJ(eonc::def::LJParams ljp)
      : Potential(PotType::LJ),
        u0{ljp.u0},
        cuttOffR{ljp.cutoff},
        psi{ljp.psi} {}

  ~LJ() = default;

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

  void setParameters(double r0Recieved, double u0Recieved, double psiRecieved);
};
