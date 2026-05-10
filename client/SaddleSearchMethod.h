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
#include "StatusTypes.h"

#include <string_view>

namespace eonc {

class SaddleSearchMethod {
protected:
  std::shared_ptr<Potential> pot;
  const Parameters &params;

public:
  SaddleSearchMethod(std::shared_ptr<Potential> potPassed,
                     const Parameters &paramsPassed)
      : pot{std::move(potPassed)},
        params{paramsPassed} {}
  virtual ~SaddleSearchMethod() = default;
  /// Run the saddle search. Returns the typed terminal status; the
  /// underlying int value is preserved for results.dat / driver
  /// round-trips via to_int(SaddleStatus).
  virtual SaddleStatus run() = 0;
  virtual double getEigenvalue() = 0;
  virtual AtomMatrix getEigenvector() = 0;
  /// Human-readable message for a status value. Default forwards to
  /// the StatusTypes.h overload; subclasses override only if they
  /// need a different vocabulary.
  virtual std::string_view describeStatus(SaddleStatus status) const {
    return statusMessage(status);
  }
  virtual SaddleStatus getStatus() const { return SaddleStatus::Good; }
  virtual int getIterationCount() const { return 0; }
  virtual int getForceCalls() const { return 0; }
};

} // namespace eonc

using eonc::SaddleSearchMethod;
