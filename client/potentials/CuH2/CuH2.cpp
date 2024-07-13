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

#include "CuH2.h"
#include <set>

namespace eonc {

void CuH2::forceImpl(const ForceInput &fip, ForceOut *efvd) {
  std::multiset<int> natmc;
  int natms[2]{0, 0}; // Always Cu, then H
  double U[1]{0};     // Fortran takes an array of size 1
  int ndim{3 * static_cast<int>(fip.nAtoms)};

  for (size_t idx{0}; idx < fip.nAtoms; ++idx) {
    natmc.insert(fip.atmnrs[idx]);
  }

#ifdef EON_CHECKS
  // Check for Copper (29) and Hydrogen (1)
  if (natmc.count(29) == 0 || natmc.count(1) == 0) {
    throw std::runtime_error("The system does not have Copper or Hydrogen, but "
                             "the CuH2 potential was requested");
  }
#endif

  // Count Copper (29) and Hydrogen (1)
  natms[0] = static_cast<int>(natmc.count(29)); // Cu
  natms[1] = static_cast<int>(natmc.count(1));  // H

#ifdef EON_CHECKS
  // Check for other atom types
  if (natms[0] + natms[1] != static_cast<int>(fip.nAtoms)) {
    throw std::runtime_error("The system has other atom types, but the CuH2 "
                             "potential was requested");
  }
#endif

  // The box only takes the diagonal (assumes cubic)
  double box_eam[]{fip.box[0], fip.box[4], fip.box[8]};

  c_force_eam(natms, ndim, box_eam, const_cast<double *>(fip.pos), efvd->F, U);
  efvd->energy = U[0];
  // *U += 697.311695; // Adjust U by a constant value, approximately minimum
  // for
  //                   // the CuH2 slab
  return;
}

} // namespace eonc
