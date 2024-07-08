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

#include "Eigen.h"
#include "Matter.h"
#include "Parameters.h"

/* Define the interface for the lowest eigenvalue determination algorithm */
class LowestEigenmode {
protected:
  // make const
  std::shared_ptr<Potential> pot;
  std::shared_ptr<Parameters> params;

public:
  // stats information
  long totalForceCalls;
  double statsTorque;
  double statsCurvature;
  double statsAngle;
  long statsRotations;
  long totalIterations; // Only set by the gpr dimer
  static const char MINMODE_DIMER[];
  static const char MINMODE_GPRDIMER[];
  static const char MINMODE_LANCZOS[];

  LowestEigenmode(std::shared_ptr<Potential> potPassed,
                  std::shared_ptr<Parameters> parameters)
      : pot{potPassed},
        params{parameters} {}
  virtual ~LowestEigenmode() {}

  // void virtual initialize(Matter const *matter, AtomMatrix displacement) = 0;
  virtual void compute(std::shared_ptr<Matter> matter,
                       AtomMatrix direction) = 0;

  virtual double getEigenvalue() = 0;
  virtual AtomMatrix getEigenvector() = 0;
};
