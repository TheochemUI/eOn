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

namespace eonc {

/// Options for Matter::compare / distance. Split out of Parameters.h so
/// Matter.h does not pull the full parameter surface.
struct StructureComparisonOptions {
  double distance_difference{0.1};
  double neighbor_cutoff{3.3};
  bool check_rotation{false};
  bool indistinguishable_atoms{true};
  double energy_difference{0.01};
  bool remove_translation{true};
};

} // namespace eonc
