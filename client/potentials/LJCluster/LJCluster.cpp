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

#include "LJCluster.h"
namespace eonc {
// LJCluster::LJCluster(double u0Recieved, double cuttOffRRecieved, double
// psiRecieved){
//     this->setParameters(u0Recieved, cuttOffRRecieved, psiRecieved);
//     return;
// }

// General Functions
void LJCluster::setParameters(double u0Recieved, double psiRecieved) {
  _lj.setParameters({.u0 = u0Recieved,
                     .cutoff_R = std::numeric_limits<double>::infinity(),
                     .psi = psiRecieved});
  return;
}

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void LJCluster::forceImpl(const ForceInput &fip, ForceOut *efvd) {
  _lj.force({.nAtoms = fip.nAtoms,
             .pos = fip.pos,
             .atmnrs = fip.atmnrs,
             .box = nullptr},
            efvd);
  return;
}

} // namespace eonc
