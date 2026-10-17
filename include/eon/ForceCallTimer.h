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
/// Usage: { ForceCallTimer timer(targetCounter, *pot); /* work */ }
/// On destruction, adds the force calls `pot` evaluated in the scope to
/// targetCounter. The potential's own counter is read, so the count holds
/// whichever registry the potential charges.
class ForceCallTimer {
  size_t &target_;
  const Potential &pot_;
  size_t initial_;

public:
  ForceCallTimer(size_t &target, const Potential &pot)
      : target_(target), pot_(pot), initial_(pot.forceCallCounter.load()) {}
  ~ForceCallTimer() { target_ += pot_.forceCallCounter.load() - initial_; }

  // Non-copyable, non-movable
  ForceCallTimer(const ForceCallTimer &) = delete;
  ForceCallTimer &operator=(const ForceCallTimer &) = delete;
};

} // namespace eonc
