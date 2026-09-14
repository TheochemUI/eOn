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

#include "eon/Matter.h"
#include "eon/Potential.h"

#include <string>

class AMS_IO : public Potential {

public:
  AMS_IO(const Parameters &p);
  ~AMS_IO();
  void initialize() {};
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);

private:
  void passToSystem(long N, const double *R, const int *atomicNrs,
                    const double *box);
  void recieveFromSystem(long N, double *F, double *U);
  // Owned by value: the Parameters these came from may go out of scope long
  // before the potential does.
  std::string engine;
  std::string model;
  std::string forcefield;
  std::string xc;
};
