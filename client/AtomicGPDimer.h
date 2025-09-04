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
#include "LowestEigenmode.h"
#include "Matter.h"
#include "Parameters.h"
#include <vector>

// Defined in the gpr_optim target
#include "subprojects/gpr_optim/gpr/AtomicDimer.h"
#include "subprojects/gpr_optim/gpr/auxiliary/ProblemSetUp.h"

// dimer method to find the lowest curvature mode
class AtomicGPDimer : public LowestEigenmode {

public:
  // Optimization for the dimer
  static const char OPT_SCG[];
  static const char OPT_LBFGS[];

  AtomicGPDimer(std::shared_ptr<Matter> matter,
                std::shared_ptr<Parameters> params,
                std::shared_ptr<Potential> pot);
  ~AtomicGPDimer() = default;

  void compute(std::shared_ptr<Matter> matter, AtomMatrix initialDirection);
  double getEigenvalue();
  AtomMatrix getEigenvector();

private:
  std::shared_ptr<Matter> matterCenter; // initial center of the dimer
  AtomMatrix direction;                 // direction along the dimer
  AtomMatrix rotationalPlane; // direction normal to the plane of dimer rotation

  gpr::InputParameters p;
  atmd::AtomicDimer atomic_dimer;
  aux::ProblemSetUp problem_setup;
  gpr::AtomsConfiguration atoms_config;
  gpr::Observation init_observations, init_middle_point;
  gpr::Coord orient_init, R_init;
};
