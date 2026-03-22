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

#include "Eigen.h"
#include "Parameters.h"
#include <variant>

namespace eonc::neb {

/// Mills, Jonsson, Schenter, Surf. Sci. 324:305, 1995.
/// Simple forward-difference tangent: uses only the direction to the next
/// image.
struct SimpleTangent {
  AtomMatrix compute(const AtomMatrix &posDiffNext,
                     const AtomMatrix &posDiffPrev, double energy,
                     double energyPrev, double energyNext) const;
};

/// Henkelman & Jonsson, JCP 113:9978, 2000.
/// Improved tangent that selects uphill direction at monotonic points and
/// uses energy-weighted interpolation at extrema.
struct ImprovedTangent {
  AtomMatrix compute(const AtomMatrix &posDiffNext,
                     const AtomMatrix &posDiffPrev, double energy,
                     double energyPrev, double energyNext) const;
};

using TangentStrategy = std::variant<SimpleTangent, ImprovedTangent>;

/// Build the tangent strategy from parameters.
TangentStrategy buildTangentStrategy(const Parameters &params);

} // namespace eonc::neb
