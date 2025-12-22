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
  XTBPot(std::shared_ptr<Parameters> p)
      : Potential(PotType::XTB, p),
        xtb_acc{p->xtb_acc},
        xtb_electronic_temperature{p->xtb_elec_temperature},
        xtb_max_iter{p->xtb_maxiter},
        total_charge{p->xtb_charge},
        uhf{p->xtb_uhf} {
    counter = 0;
    initialized = false;
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
    res = xtb_newResults();
    if (!calc) {
      xtb_delResults(&res);
      throw std::runtime_error("Failed to create xtb results");
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
    if (res) {
      xtb_delResults(&res);
    }
    if (calc) {
      xtb_delCalculator(&calc);
    }
    if (mol) {
      xtb_delMolecule(&mol);
    }
    if (env) {
      xtb_delEnvironment(&env);
    }
    SPDLOG_INFO("[XTB] called potential {} times", counter++);
  }

  // Disable copy to prevent double-free of Fortran pointers
  XTBPot(const XTBPot &) = delete;
  XTBPot &operator=(const XTBPot &) = delete;

  // To satisfy interface
  void cleanMemory(void);

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  enum class GFNMethod { GFNFF, GFN0xTB, GFN1xTB, GFN2xTB };
  xtb_TEnvironment env = nullptr;
  xtb_TCalculator calc = nullptr;
  xtb_TMolecule mol = nullptr;
  xtb_TResults res = nullptr;

  GFNMethod xtb_paramset;
  double xtb_acc;
  double xtb_electronic_temperature;
  size_t xtb_max_iter;
  double total_charge = 0.0;
  int uhf = 0;

  size_t counter;
  bool initialized;
};
