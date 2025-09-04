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

#include <iostream>
#include <math.h>
// #include "../../system_unit.h" // unit converters
#include "../../Potential.h"

/** Lennard Jones potential.*/
class LJCluster : public Potential {

private:
  //	Variables
  double u0;
  double cuttOffR;
  double psi;

  double cuttOffU;

public:
  // Functions
  // constructor
  LJCluster(std::shared_ptr<Parameters> params)
      : Potential(params),
        u0{1.0},
        cuttOffR{15.0},
        psi{1.0} {};
  // TODO: Put these in parameters
  // LJCluster(double r0Recieved, double u0Recieved, double psiRecieved);

  ~LJCluster();

  // Just to satisfy interface
  void cleanMemory();

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
  void setParameters(double r0Recieved, double u0Recieved, double psiRecieved);
};
