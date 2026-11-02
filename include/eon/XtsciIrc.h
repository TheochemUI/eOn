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

#include <vector>

namespace eonc {

/// Accepted points of one branch of the intrinsic reaction coordinate (IRC),
/// in order away from the saddle.
struct IrcBranch {
  std::vector<AtomMatrix> positions;
  std::vector<double> energies;
  std::vector<double> maxForces;
  /// rgsaddle steps taken; a step that rejects its first trial adds none
  /// of the points above.
  long steps{0};
  bool reachedMinimum{false};
};

/// Both branches of the IRC. `minus` leaves the saddle against `mode` and
/// `plus` along it, the sides process search relaxes into min1 and min2.
struct IrcPath {
  IrcBranch minus;
  IrcBranch plus;
};

/// True when this build can run followIrc: linked to rgsaddle at C ABI
/// minor 14 or newer, at compile time and at load time.
bool ircAvailable();

/// Follows the rgsaddle IRC from a converged saddle along both signs of its
/// lowest mode. The walk is in mass-weighted coordinates of the free atoms.
/// A branch stops at a point that rgsaddle reports as a minimum, or after
/// irc_max_steps steps of [Saddle Search]. rgsaddle makes that report only
/// after the first step and at a largest atomic force at or below
/// irc_converged_force.
IrcPath followIrc(const Matter &saddle, const AtomMatrix &mode,
                  const Parameters &params);

} // namespace eonc
