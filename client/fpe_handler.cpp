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
#include "fpe_handler.h"

#include <cfenv>
#include <csignal>
#include <iostream>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <float.h>
#endif

#if defined(__APPLE__) && defined(__x86_64__)
#include <xmmintrin.h>
#endif

namespace eonc {

#ifdef _WIN32
static LONG WINAPI windowsFPEHandler(EXCEPTION_POINTERS *info) {
  DWORD code = info->ExceptionRecord->ExceptionCode;
  switch (code) {
  case EXCEPTION_FLT_DIVIDE_BY_ZERO:
    std::cerr << "Floating point exception: Division by zero" << std::endl;
    break;
  case EXCEPTION_FLT_INVALID_OPERATION:
    std::cerr << "Floating point exception: Invalid operation" << std::endl;
    break;
  case EXCEPTION_FLT_OVERFLOW:
    std::cerr << "Floating point exception: Overflow" << std::endl;
    break;
  case EXCEPTION_FLT_UNDERFLOW:
    std::cerr << "Floating point exception: Underflow" << std::endl;
    break;
  case EXCEPTION_FLT_INEXACT_RESULT:
    std::cerr << "Floating point exception: Inexact result" << std::endl;
    break;
  case EXCEPTION_FLT_DENORMAL_OPERAND:
    std::cerr << "Floating point exception: Denormal operand" << std::endl;
    break;
  case EXCEPTION_FLT_STACK_CHECK:
    std::cerr << "Floating point exception: Stack check" << std::endl;
    break;
  default:
    return EXCEPTION_CONTINUE_SEARCH;
  }
  std::abort();
  return EXCEPTION_EXECUTE_HANDLER;
}
#else
static void fpe_signal_handler(int sig, siginfo_t *sip, void *scp) {
  int fe_code = sip->si_code;
  std::cerr << "Floating point exception: ";

  if (fe_code == FPE_FLTDIV)
    std::cerr << "Division by zero" << std::endl;
  else if (fe_code == FPE_FLTINV)
    std::cerr << "Invalid operation" << std::endl;
  else if (fe_code == FPE_FLTOVF)
    std::cerr << "Overflow" << std::endl;
  else
    std::cerr << "Code detected: " << fe_code << std::endl;

  std::abort();
}
#endif

void enableFPE() {
#ifdef _WIN32
  // Register Windows SEH handler for FPE reporting
  SetUnhandledExceptionFilter(windowsFPEHandler);
  // Enable floating-point exceptions on Windows
  _controlfp_s(nullptr, 0, _MCW_EM);
  _controlfp_s(nullptr, ~(_EM_ZERODIVIDE | _EM_INVALID | _EM_OVERFLOW),
               _MCW_EM);
#elif defined(__unix__)
  // Enable floating-point exceptions on Unix
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#elif defined(__APPLE__) && defined(__aarch64__)
  // Enable floating-point exceptions on ARM macOS
  fenv_t env;
  fegetenv(&env);
  env.__fpsr &= ~(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  fesetenv(&env);
#elif defined(__APPLE__) && defined(__x86_64__)
  // Enable floating-point exceptions on Intel macOS
  _MM_SET_EXCEPTION_MASK(
      _MM_MASK_MASK &
      ~(_MM_MASK_INVALID | _MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW));
#else
  std::cerr << "FPE trapping not supported on this platform." << std::endl;
#endif

#ifndef _WIN32
  // Register POSIX signal handler
  struct sigaction act;
  act.sa_sigaction = fpe_signal_handler;
  sigemptyset(&act.sa_mask);
  act.sa_flags = SA_SIGINFO;
  sigaction(SIGFPE, &act, nullptr);
#endif
}

void disableFPE() {
#ifdef _WIN32
  // Mask all floating-point exceptions (restore default behavior)
  unsigned int control;
  _controlfp_s(&control, _MCW_EM, _MCW_EM);
#elif defined(__unix__)
  fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#elif defined(__APPLE__)
  fenv_t env;
  fegetenv(&env);
#if defined(__aarch64__)
  env.__fpsr |= (FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
  fesetenv(&env);
#endif
}

void FPEHandler::eat_fpe() {
  std::lock_guard<std::mutex> lock(mutex_);
  feholdexcept(&orig_feenv);
}

void FPEHandler::restore_fpe() {
  std::lock_guard<std::mutex> lock(mutex_);
  fesetenv(&orig_feenv);
}

} // namespace eonc
