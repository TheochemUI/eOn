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

#include <fenv.h>
#include <mutex>

namespace eonc {

// Floating Point Trapping. It is platform specific!
// This causes the program to crash on divison by zero,
// invalid operations, and overflows.
void enableFPE(void);

class FPEHandler {
public:
  // Fix for gh-184, see
  // https://github.com/numpy/numpy/issues/20504#issuecomment-985542508
  void eat_fpe();
  void restore_fpe();

private:
  fenv_t orig_feenv;
  std::mutex mutex_;
};

} // namespace eonc
