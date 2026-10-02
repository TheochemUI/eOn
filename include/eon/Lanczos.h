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
#include "EonLogger.h"

#include "LowestEigenmode.h"
#include "Matter.h"
#include "Parameters.h"

#include <memory>

namespace eonc {

// Lanczos method to find the lowest curvature mode.
// Default Krylov space = all free atoms ([Lanczos] phva_atoms = All).
// Optional mobile list (or phva_atoms) is the PHVA active set: free/fixed is
// the optimizer mask and is not redefined.
class Lanczos : public LowestEigenmode {

public:
  Lanczos(std::shared_ptr<Matter> matter, const Parameters &params,
          std::shared_ptr<Potential> pot);
  ~Lanczos() = default;
  void compute(std::shared_ptr<Matter> matter,
               AtomMatrix initialDirection) override;
  /// Same as compute(matter, dir) but with an explicit mobile atom list
  /// (intersected with free flags). Eigenvector rows outside the list are 0.
  void compute(std::shared_ptr<Matter> matter, AtomMatrix initialDirection,
               const VectorXi &mobileAtoms);
  double getEigenvalue() override;
  AtomMatrix getEigenvector() override;
  /// The probe is the first Krylov product's displaced image, for the
  /// one-sided product only.
  void evaluateWithProbe(const std::shared_ptr<Matter> &centre,
                         const AtomMatrix &direction) override;

private:
  // Displaced image evaluated together with its centre by
  // evaluateWithProbe(); the first product of compute() takes it.
  std::unique_ptr<Matter> probe_;
  AtomMatrix lowestEv;
  double lowestEw{0.0};
  eonc::log::Scoped log;
};

} // namespace eonc
