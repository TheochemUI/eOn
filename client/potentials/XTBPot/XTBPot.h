//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#pragma once

#include "../../Potential.h"
#include "xtb.h"

class XTBPot final : public Potential {

public:
  // Functions
  XTBPot(std::shared_ptr<Parameters> p)
      : Potential(PotType::XTB, p) {
    env = xtb_newEnvironment();
    xtb_setVerbosity(env, XTB_VERBOSITY_MUTED);
    if (!env) {
      throw std::runtime_error("Failed to create xtb environment");
    }

    calc = xtb_newCalculator();
    if (!calc) {
      xtb_delEnvironment(&env);
      throw std::runtime_error("Failed to create xtb calculator");
    }
  }
  virtual ~XTBPot() {
    if (calc) {
      xtb_delCalculator(&calc);
    }
    if (env) {
      xtb_delEnvironment(&env);
    }
  }

  // To satisfy interface
  void cleanMemory(void);

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  xtb_TEnvironment env;
  xtb_TCalculator calc;
};
