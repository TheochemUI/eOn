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

#include "ConFileIO.h"
#include "Eigen.h"
#include "EigenmodeStrategy.h"
#include "EonLogger.h"
#include "Matter.h"
#include <filesystem>
#include <memory>
#include <optional>
#include <readcon-core.hpp>
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

/// Build stamped ConFrames for a NEB band (same metadata as writePathCon).
/// Empty on invalid path size. Does not write to disk.
[[nodiscard]] std::vector<readcon::ConFrame> pathToConFrames(
    const std::vector<std::shared_ptr<Matter>> &path,
    const std::vector<std::shared_ptr<AtomMatrix>> &tangent,
    const std::vector<std::shared_ptr<EigenmodeStrategy>> &eigenmode_solvers,
    long numImages, bool estimateEigenvalues,
    std::optional<size_t> bandIndex = std::nullopt);

/// Write a NEB band as a multi-frame .con via readcon ConFrameBuilder::clone().
[[nodiscard]] eonc::io::IoStatus writePathCon(
    const std::vector<std::shared_ptr<Matter>> &path,
    const std::vector<std::shared_ptr<AtomMatrix>> &tangent,
    const std::vector<std::shared_ptr<EigenmodeStrategy>> &eigenmode_solvers,
    long numImages, bool estimateEigenvalues, std::string filename,
    std::optional<size_t> bandIndex = std::nullopt);

} // namespace eonc::neb
