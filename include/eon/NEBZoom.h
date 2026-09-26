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

#include "Matter.h"
#include "Parameters.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace eonc::neb::zoom {

/// Inclusive image indices on the current band. Invalid when the
/// selection is a single image and cannot anchor a spline.
struct Window {
  std::size_t lo{0};
  std::size_t hi{0};
  bool valid{false};
};

/// Auto: contiguous images around the climbing image whose energy is
/// above E_ref + alpha * barrier. Manual: climbing image +/- offset.
/// A one-image auto window falls back to the manual offset.
Window selectWindow(const std::vector<double> &energy,
                    std::size_t climbingImage,
                    const neb_options_t::zoom_options_t &cfg);

/// Place every band image on the window by equal arc length.
/// Window endpoints become the new fixed band endpoints.
/// Cubic uses resamplePath; linear uses segment interpolation.
/// Returns false when the window cannot be resampled.
bool redistributePath(std::vector<std::shared_ptr<Matter>> &path, Window window,
                      neb_options_t::zoom_options_t::Interpolation how);

} // namespace eonc::neb::zoom
