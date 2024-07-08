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

// dimer method to find the lowest curvature mode
class ImprovedDimer : public LowestEigenmode {

private:
  std::shared_ptr<spdlog::logger> log;

public:
  // Optimization for the dimer
  //    static const std::string OPT_SD;
  //    static const std::string OPT_CG;
  //    static const std::string OPT_LBFGS;
  static const char OPT_SD[];
  static const char OPT_CG[];
  static const char OPT_LBFGS[];

  ImprovedDimer(std::shared_ptr<Matter> matter,
                std::shared_ptr<Parameters> params,
                std::shared_ptr<Potential> pot);
  ~ImprovedDimer() = default;

  void compute(std::shared_ptr<Matter> matter, AtomMatrix initialDirection);
  double getEigenvalue();
  AtomMatrix getEigenvector();

  std::shared_ptr<Matter> x0; // Center image
  std::shared_ptr<Matter> x1; // Forward image
  VectorType tau;             // Dimer direction
  VectorType theta;           // Dimer rotation direction
  VectorType F_R;             // Dimer rotational force
  double C_tau;               // Curvature along tau

  // parameters used for conjugate gradients
  VectorType F_R_Old;
  VectorType thetaOld;
  double a, b, gamma;
  bool init_cg;

  // variables for LBFGS
  std::vector<VectorType> s, y;
  std::vector<double> rho;
  bool init_lbfgs;
  VectorType rPrev;

  std::vector<VectorType> gradients;
  std::vector<VectorType> positions;
};
