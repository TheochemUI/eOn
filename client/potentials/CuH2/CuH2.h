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

// natms(2), ndim, U(1), R(ndim), F(ndim), box(3)
extern "C" void c_force_eam(int *natms, int ndim, double *box, double *R,
                            double *F, double *U);

class CuH2 : public Potential {

public:
  // Functions
  CuH2(std::shared_ptr<Parameters> p)
      : Potential(PotType::CUH2, p) {}

  // To satisfy interface
  void cleanMemory(void);

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
};
