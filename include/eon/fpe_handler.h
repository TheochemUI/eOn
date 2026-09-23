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

#include <cstdint>
#include <fenv.h>
#include <mutex>

namespace eonc {

// Windows x64 EXCEPTION_CONTINUE_EXECUTION reloads this MXCSR from
// CONTEXT.MxCsr. Bits 7-12 are the exception masks; bits 0-5 are sticky
// flags. _controlfp_s updates the live register only.
constexpr std::uint32_t kWindowsMxcsrExceptionMasks = 0x1F80u;
constexpr std::uint32_t maskWindowsMxcsrForContinue(std::uint32_t mxcsr) {
  return (mxcsr & ~std::uint32_t{0x3Fu}) | kWindowsMxcsrExceptionMasks;
}

// Floating Point Trapping. It is platform specific!
// This causes the program to crash on divison by zero,
// invalid operations, and overflows.
void enableFPE(void);
void disableFPE(void);

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

// eat_fpe for one scope. The destructor restores traps on every exit,
// including when the guarded call throws.
class FPEGuard {
public:
  FPEGuard() { handler_.eat_fpe(); }
  ~FPEGuard() { handler_.restore_fpe(); }
  FPEGuard(const FPEGuard &) = delete;
  FPEGuard &operator=(const FPEGuard &) = delete;

private:
  FPEHandler handler_;
};

} // namespace eonc
