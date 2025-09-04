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
#include "../../subprojects/gpr_optim/gpr/ml/GaussianProcessRegression.h"

/** Template to use if user want to provide potential. */
class GPRPotential : public Potential {

private:
  gpr::GaussianProcessRegression *gpr_model;

public:
  // Functions
  // constructor and destructor
  GPRPotential(Parameters *p);

  void registerGPRObject(gpr::GaussianProcessRegression *_gpr_model);

  // To satisfy interface
  void initialize(void);
  void cleanMemory(void);

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);
};
