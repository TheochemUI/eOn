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
#include "Matter.h"
#include "MinModeSaddleSearch.h"
#include "SaddleSearchMethod.h"
#include <vector>

namespace eonc {

class BasinHoppingSaddleSearch : public SaddleSearchMethod {
public:
  BasinHoppingSaddleSearch(std::shared_ptr<Matter> reactant,
                           std::shared_ptr<Matter> displacement,
                           std::shared_ptr<Potential> potPassed,
                           const Parameters &parametersPassed)
      : SaddleSearchMethod(potPassed, parametersPassed),
        reactant{std::make_shared<Matter>(*reactant)}, saddle{displacement} {
    eigenvector.resize(reactant->numberOfAtoms(), 3);
    eigenvector.setZero();
  }
  ~BasinHoppingSaddleSearch() = default;

  // Minimum-image (next - prev) / 2. setPositions wraps into the cell, so a
  // raw central difference is a box-length jump when a bead crosses a face.
  [[nodiscard]] static AtomMatrix
  initialDimerDirection(const Matter &image, const AtomMatrix &prev,
                        const AtomMatrix &next) {
    return image.pbc(next - prev) / 2.0;
  }

  int run(void);
  /// Quenched Metropolis weight for energy change `de`.
  /// Divides by `kB * temperature` only for an uphill hop when both are
  /// positive. Non-positive temperature or `kB` refuses that hop.
  static double metropolisProbability(double de, double kB, double temperature);
  double getEigenvalue();
  AtomMatrix getEigenvector();
  std::string_view describeStatus(int code) const override {
    // Code 1 is the Metropolis rejection, not MinMode STATUS_INIT.
    if (code == 1) {
      return "Basin hop rejected";
    }
    return MinModeSaddleSearch::statusMessage(code);
  }
  int getStatus() const override { return status; }

  double eigenvalue{0.0};
  AtomMatrix eigenvector;

  std::shared_ptr<Matter> reactant;
  std::shared_ptr<Matter> saddle;
  std::shared_ptr<Matter> product;

  int status{0};

private:
  eonc::log::Scoped log;
};

} // namespace eonc
