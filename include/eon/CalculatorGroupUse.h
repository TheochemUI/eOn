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

#include <algorithm>
#include <cstddef>
#include <format>
#include <numeric>
#include <string>
#include <vector>

namespace eonc {

/// How a job used its calculator groups: the seconds each group spent
/// inside engine calls, the systems it evaluated, and the driver's wall
/// time inside grouped requests. The ratios follow the POP parallel
/// efficiency model: parallel efficiency is load balance times
/// communication efficiency.
struct CalculatorGroupUse {
  std::vector<double> busy;    ///< seconds in engine calls, per group
  std::vector<double> systems; ///< systems evaluated, per group
  double wall{0.0};            ///< driver seconds inside grouped requests
  double singleWall{0.0};      ///< of which single-system requests
  long batches{0};             ///< grouped batch requests
  long singles{0};             ///< single-system requests (group 0 only)

  [[nodiscard]] std::size_t groups() const noexcept { return busy.size(); }

  [[nodiscard]] double meanBusy() const noexcept {
    if (busy.empty()) {
      return 0.0;
    }
    return std::accumulate(busy.begin(), busy.end(), 0.0) /
           static_cast<double>(busy.size());
  }

  [[nodiscard]] double maxBusy() const noexcept {
    return busy.empty() ? 0.0 : *std::max_element(busy.begin(), busy.end());
  }

  /// Mean over max busy time: 1 when every group carried the same work.
  [[nodiscard]] double loadBalance() const noexcept {
    const double m = maxBusy();
    return m > 0.0 ? meanBusy() / m : 1.0;
  }

  /// Max busy over wall: the share of the driver's wait that the busiest
  /// group spent computing. Waits for a slower group in one batch, serial
  /// requests on group 0 and the result broadcasts all lower it.
  [[nodiscard]] double communicationEfficiency() const noexcept {
    return wall > 0.0 ? std::min(1.0, maxBusy() / wall) : 1.0;
  }

  [[nodiscard]] double parallelEfficiency() const noexcept {
    return loadBalance() * communicationEfficiency();
  }

  /// Fraction of the grouped wall time a group spent outside engine calls.
  [[nodiscard]] double idleFraction(std::size_t group) const noexcept {
    if (group >= busy.size() || wall <= 0.0) {
      return 0.0;
    }
    return std::clamp(1.0 - busy[group] / wall, 0.0, 1.0);
  }

  [[nodiscard]] std::string table() const {
    std::string out = std::format(
        "calculator groups: {} groups, {} batches and {} single requests, "
        "{:.1f} s grouped wall ({:.1f} s single)\n",
        groups(), batches, singles, wall, singleWall);
    out += std::format("{:>6} {:>8} {:>10} {:>6}\n", "group", "systems",
                       "busy_s", "idle");
    for (std::size_t g = 0; g < busy.size(); ++g) {
      const double n = g < systems.size() ? systems[g] : 0.0;
      out += std::format("{:>6} {:>8.0f} {:>10.2f} {:>5.1f}%\n", g, n, busy[g],
                         100.0 * idleFraction(g));
    }
    out += std::format("load balance {:.3f}, communication efficiency {:.3f}, "
                       "parallel efficiency {:.3f}\n",
                       loadBalance(), communicationEfficiency(),
                       parallelEfficiency());
    return out;
  }
};

} // namespace eonc
