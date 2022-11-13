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

CuH2::CuH2(Parameters *p) {}

void CuH2::initialize(void) { return; }

void CuH2::cleanMemory(void) { return; }

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// address to supercell size
void CuH2::force(long N, const double *R, const int *atomicNrs, double *F,
                 double *U, const double *box, int nImages = 1) {
  std::multiset<double> natmc;
  int natms[2]{0, 0}; // Always Cu, then H
  int ndim{3 * static_cast<int>(N)};
  for (auto idx{0}; idx < N; idx++) {
    natmc.insert(atomicNrs[idx]);
  }
  natms[0] = natmc.count(29); // Cu
  natms[1] = natmc.count(1);  // H

  // The box only takes the diagonal (assumes cubic)
  double box_eam[]{box[0], box[4], box[8]};

  c_force_eam(natms, ndim, box_eam, const_cast<double *>(R), F, U);

  // for(int i=0; i<N; i++){
  //     std::cout<<F[ 3*i ]<<" "<<F[3*i+1]<<" "<<F[3*i+2]<<"\n";
  // }
  return;
}
