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
class GPRPotential final : public Potential {

private:
  std::shared_ptr<gpr::GaussianProcessRegression> m_gprm;

public:
  // Functions
  // constructor and destructor
  GPRPotential(std::shared_ptr<Parameters> a_params)
      : Potential(a_params),
        m_gprm{nullptr} {};

  void
  registerGPRObject(std::shared_ptr<gpr::GaussianProcessRegression> a_gprm);

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
};
#endif
