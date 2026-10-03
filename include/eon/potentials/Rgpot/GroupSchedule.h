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
#include <cstdint>
#include <limits>
#include <map>
#include <vector>

namespace eonc {

/// Which calculator group evaluates which system.
///
/// A group keeps the converged orbitals of every system it evaluated, under
/// the system's key (a NEB image, a ring bead). Sending a key back to the
/// group that holds its orbitals is what makes the next SCF warm, so a key
/// stays on its home group. A batch caps each group at ceil(M / G) systems:
/// a key whose home group is full moves to the least-loaded group, which
/// becomes its home, because one extra evaluation on a full group costs a
/// whole SCF while a move costs the difference between a warm and a
/// cooler start. A single request goes to the group that last evaluated
/// the geometry nearest to it.
///
/// Every decision is made on the driver and sent to the other ranks, so
/// the state here lives on the driver only.
class GroupSchedule {
public:
  explicit GroupSchedule(int groups) : groups_{std::max(groups, 1)} {}

  [[nodiscard]] int groups() const noexcept { return groups_; }

  /// The home of a key never seen: key mod G for a key of zero or more,
  /// the batch position mod G otherwise.
  [[nodiscard]] int defaultHome(std::int64_t key, long position) const {
    const std::int64_t base = key >= 0 ? key : position;
    return static_cast<int>(((base % groups_) + groups_) % groups_);
  }

  /// Groups for a batch of systems with these keys.
  [[nodiscard]] std::vector<int>
  assign(const std::vector<std::int64_t> &keys) {
    const long m = static_cast<long>(keys.size());
    const long cap = (m + groups_ - 1) / groups_;
    std::vector<long> load(static_cast<size_t>(groups_), 0);
    std::vector<int> route(static_cast<size_t>(m), -1);
    for (long j = 0; j < m; ++j) {
      const int home = homeOf(keys[static_cast<size_t>(j)], j);
      if (load[static_cast<size_t>(home)] < cap) {
        route[static_cast<size_t>(j)] = home;
        ++load[static_cast<size_t>(home)];
      }
    }
    for (long j = 0; j < m; ++j) {
      if (route[static_cast<size_t>(j)] >= 0)
        continue;
      const auto least = std::min_element(load.begin(), load.end());
      const int g = static_cast<int>(least - load.begin());
      route[static_cast<size_t>(j)] = g;
      ++*least;
    }
    for (long j = 0; j < m; ++j)
      home_[keys[static_cast<size_t>(j)]] = route[static_cast<size_t>(j)];
    return route;
  }

  /// Remember the geometry a key was last evaluated at.
  void record(std::int64_t key, const double *positions, long n3) {
    auto &stored = last_[key];
    stored.assign(positions, positions + n3);
  }

  /// The group for a single request: the home of the recorded geometry
  /// nearest to `positions` (root-mean-square over coordinates), group 0
  /// when nothing of this size is recorded.
  [[nodiscard]] int nearest(const double *positions, long n3) const {
    double best = std::numeric_limits<double>::infinity();
    int group = 0;
    for (const auto &[key, stored] : last_) {
      if (static_cast<long>(stored.size()) != n3)
        continue;
      double d2 = 0.0;
      for (long i = 0; i < n3; ++i) {
        const double d = stored[static_cast<size_t>(i)] - positions[i];
        d2 += d * d;
      }
      const auto home = home_.find(key);
      if (d2 < best && home != home_.end()) {
        best = d2;
        group = home->second;
      }
    }
    return group;
  }

  /// A single request's group becomes the home of the single-request key.
  void settle(std::int64_t key, int group) { home_[key] = group; }

private:
  [[nodiscard]] int homeOf(std::int64_t key, long position) const {
    const auto it = home_.find(key);
    return it != home_.end() ? it->second : defaultHome(key, position);
  }

  int groups_;
  std::map<std::int64_t, int> home_;
  std::map<std::int64_t, std::vector<double>> last_;
};

} // namespace eonc
