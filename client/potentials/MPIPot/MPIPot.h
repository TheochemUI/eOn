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
#include "../../Parameters.h"
#include "../../Potential.h"
namespace eonc {
class MPIPot : public Potential {

public:
  MPIPot(Parameters &p);
  ~MPIPot();
  void initialize() {};
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);

private:
  int potentialRank;
  double poll_period;
};
} // namespace eonc
