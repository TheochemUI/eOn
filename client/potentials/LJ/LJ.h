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
#include <cmath>

/// Lennard-Jones potential.
class LJ final : public Potential {
private:
  double u0;
  double cuttOffR{0.0};
  double psi{0.0};
  double cuttOffU{0.0};

public:
  explicit LJ(const Parameters &params)
      : LJ(PotType::LJ, params, 1.0, 15.0, 1.0) {}

  LJ(PotType ptype, const Parameters &params, double u0_in, double cutoff_in,
     double psi_in)
      : Potential(ptype, params), u0{u0_in}, cuttOffR{cutoff_in}, psi{psi_in} {
    double r6 = std::pow(psi / cuttOffR, 6);
    cuttOffU = 4.0 * u0 * r6 * (r6 - 1.0);
  }

  ~LJ() override = default;

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

  void setParameters(double u0In, double cutoffIn, double psiIn);
};
