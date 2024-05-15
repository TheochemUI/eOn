//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "CuH2.h"
#include <set>

void CuH2::cleanMemory(void) { return; }

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// address to supercell size
void CuH2::force(long N, const double *R, const int *atomicNrs, double *F,
                 double *U, double *variance, const double *box) {
  variance = nullptr;
  std::multiset<int> natmc;
  int natms[2]{0, 0}; // Always Cu, then H
  int ndim{3 * static_cast<int>(N)};

  for (long idx{0}; idx < N; ++idx) {
    natmc.insert(atomicNrs[idx]);
  }

#ifndef NDEBUG
  // Check for Copper (29) and Hydrogen (1)
  if (natmc.count(29) == 0 || natmc.count(1) == 0) {
    throw std::runtime_error("The system does not have Copper or Hydrogen, but "
                             "the CuH2 potential was requested");
  }
#endif

  // Count Copper (29) and Hydrogen (1)
  natms[0] = static_cast<int>(natmc.count(29)); // Cu
  natms[1] = static_cast<int>(natmc.count(1));  // H

#ifndef NDEBUG
  // Check for other atom types
  if (natms[0] + natms[1] != N) {
    throw std::runtime_error("The system has other atom types, but the CuH2 "
                             "potential was requested");
  }
#endif

  // The box only takes the diagonal (assumes cubic)
  double box_eam[]{box[0], box[4], box[8]};

  c_force_eam(natms, ndim, box_eam, const_cast<double *>(R), F, U);
  *U += 697.311695; // Adjust U by a constant value, approximately minimum for
                    // the CuH2 slab
  return;
}
