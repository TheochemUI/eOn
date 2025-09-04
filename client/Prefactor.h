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

namespace Prefactor {
const char RATE_HTST[] = "htst";
const char RATE_QQHTST[] = "qqhtst";
const char FILTER_CUTOFF[] = "cutoff";
const char FILTER_FRACTION[] = "fraction";

int getPrefactors(Parameters *parameters, Matter *min1, Matter *saddle,
                  Matter *min2, double &pref1, double &pref2);
VectorXi movedAtoms(Parameters *parameters, Matter *min1, Matter *saddle,
                    Matter *min2);
VectorXi movedAtomsPct(Parameters *parameters, Matter *min1, Matter *saddle,
                       Matter *min2);
VectorXi allFreeAtoms(Matter *matter);
VectorXd removeZeroFreqs(Parameters *parameters, VectorXd freqs);
void logFreqs(VectorXd freqs, char *name);
} // namespace Prefactor
