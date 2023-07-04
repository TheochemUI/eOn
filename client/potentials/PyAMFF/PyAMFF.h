//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef PYAMFF_POTENTIAL
#define PYAMFF_POTENTIAL

#include "../../Potential.h"

// extern "C" void calc_eon(long *nAtoms, const double [], const double [],
// const int [], double *F, double *U); extern "C" {
//     void calceon_(const long int *nAtoms, const double *R, const double *box,
//     const int *atomicNrs, double *F, double *U);
// }

class PyAMFF : public Potential {

public:
  PyAMFF(void);
  ~PyAMFF();
  void initialize(){};
  void cleanMemory(void);
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);
  bool new_pyamff;
};
#endif
