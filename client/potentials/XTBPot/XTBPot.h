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
#include "units.hpp"
#include "xtb.h"

class XTBPot final : public Potential {
public:
  // Functions
  XTBPot(std::shared_ptr<Parameters> p)
      : Potential(PotType::XTB, p),
        xtb_acc{p->xtb_acc},
        xtb_electronic_temperature{p->xtb_elec_temperature},
        xtb_max_iter{p->xtb_maxiter} {
    counter = 0;
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
    // Unmarshal parameters
    if (p->xtb_paramset == "GFNFF") {
      xtb_paramset = GFNMethod::GFNFF;
    } else if (p->xtb_paramset == "GFN0xTB") {
      xtb_paramset = GFNMethod::GFN0xTB;
    } else if (p->xtb_paramset == "GFN1xTB") {
      xtb_paramset = GFNMethod::GFN1xTB;
    } else if (p->xtb_paramset == "GFN2xTB") {
      xtb_paramset = GFNMethod::GFN2xTB;
    } else {
      throw std::runtime_error("Parameter set for XTB must be one of GFNFF, "
                               "GFN0xTB, GFN1xTB or GFN2xTB.\n");
    }
  }

  virtual ~XTBPot() {
    if (calc) {
      xtb_delCalculator(&calc);
    }
    if (env) {
      xtb_delEnvironment(&env);
    }
    SPDLOG_INFO("[XTB] called potential {} times", counter++);
  }

  // To satisfy interface
  void cleanMemory(void);

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  enum class GFNMethod { GFNFF, GFN0xTB, GFN1xTB, GFN2xTB };
  xtb_TEnvironment env;
  xtb_TCalculator calc;
  GFNMethod xtb_paramset;
  double xtb_acc;
  double xtb_electronic_temperature;
  size_t xtb_max_iter;
  size_t counter;
};
