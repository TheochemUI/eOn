//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef CUH2_INTERFACE
#define CUH2_INTERFACE

#include "../../Potential.h"

// natms(2), ndim, U(1), R(ndim), F(ndim), box(3)
extern "C" void c_force_eam(int *natms, int ndim, double *box, double *R,
                            double *F, double *U);

class CuH2 : public Potential {

private:
public:
  // Functions
  // constructor and destructor
  CuH2(Parameters *p);

  // To satisfy interface
  void cleanMemory(void);

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, const double *box);
  std::pair<double, AtomMatrix> get_ef(const AtomMatrix pos,
                                       const VectorXi atmnrs,
                                       const Matrix3d box) override {
    double energy{std::numeric_limits<double>::infinity()};
    long nAtoms{pos.rows()};
    AtomMatrix forces{Eigen::MatrixXd::Zero(nAtoms, 3)};
    this->force(nAtoms, pos.data(), atmnrs.data(), forces.data(), &energy,
                box.data());
    return std::make_pair(energy, forces);
  };
};
#endif
