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
#include "Potential.h"

namespace eonc {

/// RAII wrapper for tracking force calls over a scope.
/// Usage: { ForceCallTimer timer(targetCounter); /* work */ }
/// On destruction, adds the delta force calls to targetCounter.
/// The delta is the change in the sum of every live registry.
class ForceCallTimer {
  size_t &target_;
  size_t initial_;

public:
  explicit ForceCallTimer(size_t &target)
      : target_(target), initial_(PotRegistry::sum_force_calls()) {}
  ~ForceCallTimer() {
    const size_t now = PotRegistry::sum_force_calls();
    if (now >= initial_) {
      target_ += now - initial_;
    }
  }

  // Non-copyable, non-movable
  ForceCallTimer(const ForceCallTimer &) = delete;
  ForceCallTimer &operator=(const ForceCallTimer &) = delete;
};

} // namespace eonc
