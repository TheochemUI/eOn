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

#include <string_view>

namespace eonc {

class SaddleSearchMethod {
protected:
  std::shared_ptr<Potential> pot;
  const Parameters &params;

public:
  SaddleSearchMethod(std::shared_ptr<Potential> potPassed,
                     const Parameters &paramsPassed)
      : pot{potPassed},
        params{paramsPassed} {};
  virtual ~SaddleSearchMethod() {};
  virtual int run() = 0;
  virtual double getEigenvalue() = 0;
  virtual AtomMatrix getEigenvector() = 0;
  virtual std::string_view describeStatus(int status) const = 0;

  int status{0};
  int iteration{0};  // Number of iterations (for reporting)
  int forcecalls{0}; // Force calls during saddle search
};

} // namespace eonc

using eonc::SaddleSearchMethod;
