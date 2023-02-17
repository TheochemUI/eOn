//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include"Tersoff.h"

// pointer to number of atoms, pointer to array of positions
// pointer to array of forces, pointer to internal energy
// address to supercell size
std::pair<double, AtomMatrix> Tersoff::get_ef(const AtomMatrix pos,
                                             const VectorXi atmnrs,
                                             const Matrix3d m_box) {
  double energy{0};
  long N{pos.rows()};
  AtomMatrix forces{Eigen::MatrixXd::Zero(N, 3)};
  const double *R = pos.data();
  const double *box = m_box.data();
  const int *atomicNrs = atmnrs.data();
  double *F = forces.data();
  double *U = &energy;

  tersoff_(&N, R, F, U, &box[0], &box[4], &box[8]);
  return std::make_pair(energy, forces);
}
