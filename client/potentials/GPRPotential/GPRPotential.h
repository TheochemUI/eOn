//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef GPRPOT_INTERFACE
#define GPRPOT_INTERFACE

#include "../../Potential.h"
#include "../../subprojects/gpr_optim/gpr/ml/GaussianProcessRegression.h"

/** Template to use if user want to provide potential. */
class GPRPotential : public Potential {

private:
  gpr::GaussianProcessRegression *gpr_model;

public:
  // Functions
  // constructor and destructor
  GPRPotential(Parameters *p);

  void registerGPRObject(gpr::GaussianProcessRegression *_gpr_model);

  // To satisfy interface
  void initialize(void);
  void cleanMemory(void);

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box);
};
#endif
