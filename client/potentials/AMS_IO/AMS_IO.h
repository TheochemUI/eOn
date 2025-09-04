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

#include "../../Matter.h"
#include "../../Potential.h"

class AMS_IO : public Potential {

public:
  AMS_IO(std::shared_ptr<Parameters> p);
  ~AMS_IO();
  void initialize() {};
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);

private:
  void passToSystem(long N, const double *R, const int *atomicNrs,
                    const double *box);
  void recieveFromSystem(long N, double *F, double *U);
  const char *engine;
  const char *model;
  const char *forcefield;
  const char *xc;
};
