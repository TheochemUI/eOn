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
namespace eonc {
/** Lennard Jones potential.*/
class LJCluster : public Potential {

private:
  // Variables
  double u0;
  double cuttOffR;
  double psi;

  double cuttOffU;

public:
  // Functions
  // constructor
  LJCluster(eonc::def::LJParams &ljp)
      : Potential(PotType::LJCLUSTER),
        u0{ljp.u0},
        cuttOffR{ljp.cutoff},
        psi{ljp.psi} {}

  ~LJCluster();

  // Just to satisfy interface
  void cleanMemory();

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
  void setParameters(double r0Recieved, double u0Recieved, double psiRecieved);
};
} // namespace eonc
