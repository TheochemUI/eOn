/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/

#pragma once

#include "../../Potential.h"
#include "units.hpp"
#include "xtb.h"

class XTBPot final : public Potential {
public:
  // Functions
  XTBPot(Parameters &a_p)
      : Potential(PotType::XTB, a_p),
        xtb_acc{a_p.xtbpot.acc},
        xtb_electronic_temperature{a_p.xtbpot.elec_temperature},
        xtb_max_iter{a_p.xtbpot.maxiter} {
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
    if (a_p.xtbpot.paramset == "GFNFF") {
      xtb_paramset = GFNMethod::GFNFF;
    } else if (a_p.xtbpot.paramset == "GFN0xTB") {
      xtb_paramset = GFNMethod::GFN0xTB;
    } else if (a_p.xtbpot.paramset == "GFN1xTB") {
      xtb_paramset = GFNMethod::GFN1xTB;
    } else if (a_p.xtbpot.paramset == "GFN2xTB") {
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
