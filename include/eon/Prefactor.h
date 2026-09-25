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
#include "Matter.h"

#include "Parameters.h"

#include <string_view>

namespace eonc {

namespace Prefactor {
inline constexpr std::string_view RATE_HTST{"htst"};
inline constexpr std::string_view RATE_QQHTST{"qqhtst"};
inline constexpr std::string_view FILTER_CUTOFF{"cutoff"};
inline constexpr std::string_view FILTER_FRACTION{"fraction"};

int getPrefactors(const Parameters &parameters, Matter *min1, Matter *saddle,
                  Matter *min2, double &pref1, double &pref2);
VectorXi movedAtoms(const Parameters &parameters, Matter *min1, Matter *saddle,
                    Matter *min2);
VectorXi movedAtomsPct(const Parameters &parameters, Matter *min1,
                       Matter *saddle, Matter *min2);
VectorXi allFreeAtoms(Matter *matter);
void logFreqs(const VectorXd &freqs, const char *name);
} // namespace Prefactor

} // namespace eonc
