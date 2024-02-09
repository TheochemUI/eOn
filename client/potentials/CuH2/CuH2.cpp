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

void CuH2::peturb_positions(AtomMatrix &positions,
                            const Eigen::VectorXi &atmNumVec,
                            const double hcu_dist, const double hh_dist) {
  std::vector<int> hIndices, cuIndices;
  for (int i = 0; i < atmNumVec.size(); ++i) {
    if (atmNumVec[i] == 1) { // Hydrogen atom
      hIndices.push_back(i);
    } else if (atmNumVec[i] == 29) { // Copper atom
      cuIndices.push_back(i);
    } else {
      throw std::runtime_error("Unexpected atomic number");
    }
  }

  if (hIndices.size() != 2) {
    throw std::runtime_error("Expected exactly two hydrogen atoms");
  }

  // Compute the midpoint of the hydrogens
  Eigen::VectorXd hMidpoint =
      (positions.row(hIndices[0]) + positions.row(hIndices[1])) / 2;

  // Compute the HH direction
  Eigen::VectorXd hh_direction;
  if (positions(hIndices[0], 0) < positions(hIndices[1], 0)) {
    hh_direction =
        (positions.row(hIndices[1]) - positions.row(hIndices[0])).normalized();
  } else {
    hh_direction =
        (positions.row(hIndices[0]) - positions.row(hIndices[1])).normalized();
  }

  // Set the new position of the hydrogens
  positions.row(hIndices[0]) = hMidpoint - (0.5 * hh_dist) * hh_direction;
  positions.row(hIndices[1]) = hMidpoint + (0.5 * hh_dist) * hh_direction;

  // Find the z-coordinate of the topmost Cu layer
  double maxCuZ = std::numeric_limits<double>::lowest();
  for (int cuIndex : cuIndices) {
    maxCuZ = std::max(maxCuZ, positions(cuIndex, 2));
  }

  // Compute the new z-coordinate for the H atoms
  double new_z = maxCuZ + hcu_dist;

  // Update the z-coordinates of the H atoms
  for (int hIndex : hIndices) {
    positions(hIndex, 2) = new_z;
  }
}
