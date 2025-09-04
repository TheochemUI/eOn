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

// extern "C" void calc_eon(long *nAtoms, const double [], const double [],
// const int [], double *F, double *U); extern "C" {
//     void calceon_(const long int *nAtoms, const double *R, const double *box,
//     const int *atomicNrs, double *F, double *U);
// }

class PyAMFF : public Potential {

public:
  PyAMFF(void);
  ~PyAMFF();
  void initialize() {};
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);
  bool new_pyamff;
};
