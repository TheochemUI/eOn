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
class ForceCallTimer {
  size_t &target_;
  size_t initial_;

public:
  explicit ForceCallTimer(size_t &target)
      : target_(target),
        initial_(PotRegistry::get().total_force_calls()) {}
  ~ForceCallTimer() {
    target_ += PotRegistry::get().total_force_calls() - initial_;
  }

  // Non-copyable, non-movable
  ForceCallTimer(const ForceCallTimer &) = delete;
  ForceCallTimer &operator=(const ForceCallTimer &) = delete;
};

} // namespace eonc
