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

#include "Parameters.h"
#include "Potential.h"

class SaddleSearchMethod {
protected:
  std::shared_ptr<Potential> pot;
  std::shared_ptr<Parameters> params;

public:
  SaddleSearchMethod(std::shared_ptr<Potential> potPassed,
                     std::shared_ptr<Parameters> paramsPassed)
      : pot{potPassed},
        params{paramsPassed} {};
  virtual ~SaddleSearchMethod() {};
  virtual int run() = 0;
  virtual double getEigenvalue() = 0;
  virtual AtomMatrix getEigenvector() = 0;

  int status;
};
