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
#define PREFACTOR_H

#include "Eigen.h"
#include "Matter.h"

#include "Parameters.h"
namespace eonc {
namespace Prefactor {
int getPrefactors(Parameters *parameters, Matter *min1, Matter *saddle,
                  Matter *min2, double &pref1, double &pref2);
Vector<int> movedAtoms(Parameters *parameters, Matter *min1, Matter *saddle,
                       Matter *min2);
Vector<int> movedAtomsPct(Parameters *parameters, Matter *min1, Matter *saddle,
                          Matter *min2);
Vector<int> allFreeAtoms(Matter *matter);
VectorType removeZeroFreqs(Parameters *parameters, VectorType freqs);
void logFreqs(VectorType freqs, char *name);
} // namespace Prefactor

} // namespace eonc
