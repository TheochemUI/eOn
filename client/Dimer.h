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
#include "LowestEigenmode.h"
#include "Matter.h"
#include "Parameters.h"
#include "Potential.h"

namespace eonc {

/// Classic dimer method to find the lowest curvature mode.
/// Uses finite-difference rotation to converge on the minimum eigenmode.
class Dimer : public LowestEigenmode {
public:
  Dimer(const std::shared_ptr<Matter> &matter, const Parameters &params,
        const std::shared_ptr<Potential> &pot);
  ~Dimer() = default;

  void compute(const std::shared_ptr<Matter> &matter,
               AtomMatrix initialDirection);
  [[nodiscard]] double getEigenvalue();
  [[nodiscard]] AtomMatrix getEigenvector();

private:
  eonc::log::FileScoped log{"dimer", "dimer.log"};
  std::shared_ptr<Matter> matterCenter;
  std::shared_ptr<Matter> matterDimer;
  AtomMatrix direction;
  AtomMatrix rotationalPlane;
  double eigenvalue{0.0};
  long nAtoms{0};

  /// Determine rotational plane via conjugate gradient.
  void determineRotationalPlane(const AtomMatrix &rotationalForce,
                                AtomMatrix &rotationalForceOld,
                                const AtomMatrix &rotationalPlaneOld,
                                double &lengthRotationalForceOld);

  /// Rotate the dimer by the given angle (radians).
  void rotate(double rotationAngle);

  /// Compute rotational force and return curvature along the dimer.
  double calcRotationalForceReturnCurvature(AtomMatrix &rotationalForce);
};

} // namespace eonc

using eonc::Dimer;
