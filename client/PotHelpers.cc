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
#include "PotHelpers.hpp"

namespace eonc::pot {

// Typically this is done by the caller, however it is here as a sanity check
void zeroForceOut(const size_t &nAtoms, ForceOut *efvd) {
  efvd->energy = 0;
  efvd->variance = 0;
  for (size_t idx{0}; idx < nAtoms; idx++) {
    efvd->F[3 * idx] = 0;
    efvd->F[3 * idx + 1] = 0;
    efvd->F[3 * idx + 2] = 0;
  }
};

} // namespace eonc::pot
