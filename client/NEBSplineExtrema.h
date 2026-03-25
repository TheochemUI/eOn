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
#include "EigenmodeStrategy.h"
#include "EonLogger.h"
#include "Matter.h"
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace eonc::neb {

struct ExtremaResult {
  long numExtrema{0};
  std::vector<double> positions;
  std::vector<double> energies;
  std::vector<double> curvatures;
};

/// Find extrema along the MEP using cubic spline interpolation.
ExtremaResult
findSplineExtrema(const std::vector<std::shared_ptr<Matter>> &path,
                  const std::vector<std::shared_ptr<AtomMatrix>> &tangent,
                  long numImages);

/// Print NEB image data to log and optionally to file.
void printImageData(
    const std::vector<std::shared_ptr<Matter>> &path,
    const std::vector<std::shared_ptr<AtomMatrix>> &tangent,
    const std::vector<std::shared_ptr<EigenmodeStrategy>> &eigenmode_solvers,
    long numImages, bool estimateEigenvalues, bool writeToFile, size_t idx,
    eonc::log::Scoped log);

} // namespace eonc::neb
