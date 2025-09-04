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

/** Template to use if user want to provide potential. */
class NewPot : public Potential {

private:
  //	Variables
  double fake1;
  double fake2;

public:
  // Functions
  // constructor and destructor
  NewPot(std::shared_ptr<Parameters> p)
      : Potential(p),
        fake1{0},
        fake2{0} {};

  // To satisfy interface
  void cleanMemory(void);

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
};
