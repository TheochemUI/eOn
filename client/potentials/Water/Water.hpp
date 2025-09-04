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
/** @file
Wrapper for Eon
@author Jean-Claude C. Berthet
@date 2007
University of Iceland
*/

#include "../../Potential.h"
#include "spce_ccl.hpp"
#include "tip4p_ccl.hpp"

class Tip4p : public Potential, private forcefields::Tip4p {
public:
  Tip4p(std::shared_ptr<Parameters> params)
      : Potential(params),
        forcefields::Tip4p(8.5, 1.0) {};
  // Functions
  // constructor and destructor

  // To satisfy interface
  void cleanMemory(void) {}
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, const double *box) override;
};

class SpceCcl : public Potential, private forcefields::SpceCcl {
public:
  SpceCcl(std::shared_ptr<Parameters> params)
      : Potential(params),
        forcefields::SpceCcl(8.5, 1.0) {}
  // Functions
  // constructor and destructor

  // To satisfy interface
  void cleanMemory(void) {}
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;
};
