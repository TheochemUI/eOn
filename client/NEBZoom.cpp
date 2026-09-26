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
#include "eon/NEBZoom.h"
#include "eon/NEBInitialPaths.hpp"

#include <algorithm>
#include <cmath>
#include <span>
#include <utility>

namespace eonc::neb::zoom {
namespace {

std::pair<std::size_t, std::size_t>
manualWindow(std::size_t n, std::size_t climbingImage, int offset) {
  const auto half = static_cast<std::size_t>(std::max(offset, 1));
  const std::size_t lo = climbingImage > half ? climbingImage - half : 0;
  const std::size_t hi = std::min(n - 1, climbingImage + half);
  return {lo, hi};
}

std::vector<Matter> linearResample(const std::vector<Matter> &sub,
                                   std::size_t count) {
  const std::size_t n = sub.size();
  std::vector<double> arc(n, 0.0);
  for (std::size_t i = 1; i < n; ++i) {
    AtomMatrix diff =
        sub[i].pbc(sub[i].getPositions() - sub[i - 1].getPositions());
    arc[i] = arc[i - 1] + diff.norm();
  }
  const double total = arc.back();
  std::vector<Matter> placed;
  placed.reserve(count);
  if (count == 0) {
    return placed;
  }
  placed.push_back(sub.front());
  for (std::size_t i = 1; i + 1 < count; ++i) {
    const double target = (total > 1e-12) ? total * static_cast<double>(i) /
                                                static_cast<double>(count - 1)
                                          : 0.0;
    std::size_t lo = 0;
    for (std::size_t j = 1; j < n; ++j) {
      if (arc[j] >= target) {
        lo = j - 1;
        break;
      }
    }
    const std::size_t hi = std::min(lo + 1, n - 1);
    const double span = arc[hi] - arc[lo];
    const double f = (span > 1e-12) ? (target - arc[lo]) / span : 0.0;
    Matter image(sub.front());
    image.setPositions((1.0 - f) * sub[lo].getPositions() +
                       f * sub[hi].getPositions());
    placed.push_back(std::move(image));
  }
  if (count > 1) {
    placed.push_back(sub.back());
  }
  return placed;
}

} // namespace

Window selectWindow(const std::vector<double> &energy,
                    std::size_t climbingImage,
                    const neb_options_t::zoom_options_t &cfg) {
  Window window;
  const std::size_t n = energy.size();
  if (n < 3 || climbingImage == 0 || climbingImage + 1 >= n) {
    return window;
  }

  std::size_t lo = 0;
  std::size_t hi = 0;
  if (cfg.mode == neb_options_t::zoom_options_t::Mode::Manual) {
    std::tie(lo, hi) = manualWindow(n, climbingImage, cfg.offset);
  } else {
    const double eRef = std::min(energy.front(), energy.back());
    const double eMax = *std::max_element(energy.begin(), energy.end());
    const double barrier = eMax - eRef;
    const bool usable = barrier > 0.0 && cfg.alpha > 0.0 && cfg.alpha < 1.0;
    if (!usable) {
      std::tie(lo, hi) = manualWindow(n, climbingImage, cfg.offset);
    } else {
      const double threshold = eRef + cfg.alpha * barrier;
      lo = climbingImage;
      hi = climbingImage;
      while (lo > 0 && energy[lo - 1] > threshold) {
        --lo;
      }
      while (hi + 1 < n && energy[hi + 1] > threshold) {
        ++hi;
      }
      if (hi <= lo) {
        std::tie(lo, hi) = manualWindow(n, climbingImage, cfg.offset);
      }
    }
  }
  if (hi <= lo || hi >= n) {
    return window;
  }
  window.lo = lo;
  window.hi = hi;
  window.valid = true;
  return window;
}

bool redistributePath(std::vector<std::shared_ptr<Matter>> &path, Window window,
                      neb_options_t::zoom_options_t::Interpolation how) {
  if (!window.valid || path.size() < 3 || window.hi <= window.lo ||
      window.hi >= path.size()) {
    return false;
  }
  for (const auto &image : path) {
    if (!image) {
      return false;
    }
  }

  std::vector<Matter> sub;
  sub.reserve(window.hi - window.lo + 1);
  for (std::size_t i = window.lo; i <= window.hi; ++i) {
    sub.push_back(*path[i]);
  }

  std::vector<Matter> placed;
  if (how == neb_options_t::zoom_options_t::Interpolation::Linear) {
    placed = linearResample(sub, path.size());
  } else if (sub.size() == path.size()) {
    std::vector<std::shared_ptr<Matter>> tmp;
    tmp.reserve(sub.size());
    for (const auto &image : sub) {
      tmp.push_back(std::make_shared<Matter>(image));
    }
    eonc::helpers::neb_paths::resamplePathInPlace(
        std::span<std::shared_ptr<Matter>>{tmp.data(), tmp.size()});
    placed.reserve(tmp.size());
    for (const auto &image : tmp) {
      placed.push_back(*image);
    }
  } else {
    placed = eonc::helpers::neb_paths::resamplePath(sub, path.size() - 2);
  }
  if (placed.size() != path.size()) {
    return false;
  }
  for (std::size_t i = 0; i < path.size(); ++i) {
    path[i]->setPositions(placed[i].getPositions());
  }
  return true;
}

} // namespace eonc::neb::zoom
