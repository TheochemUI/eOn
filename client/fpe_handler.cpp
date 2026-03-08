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
#include <cstdio>

#ifndef _WIN32
#include <unistd.h>
#endif

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <float.h>
#include <windows.h>
#endif

#if defined(__linux__)
#include <ucontext.h>
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
    fprintf(stderr, "FPE (continuing): division by zero\n");
    break;
  case EXCEPTION_FLT_INVALID_OPERATION:
    fprintf(stderr, "FPE (continuing): invalid operation\n");
    break;
  case EXCEPTION_FLT_OVERFLOW:
    fprintf(stderr, "FPE (continuing): overflow\n");
    break;
  case EXCEPTION_FLT_UNDERFLOW:
    fprintf(stderr, "FPE (continuing): underflow\n");
    break;
  case EXCEPTION_FLT_INEXACT_RESULT:
    fprintf(stderr, "FPE (continuing): inexact result\n");
    break;
  case EXCEPTION_FLT_DENORMAL_OPERAND:
    fprintf(stderr, "FPE (continuing): denormal operand\n");
    break;
  case EXCEPTION_FLT_STACK_CHECK:
    fprintf(stderr, "FPE (continuing): stack check\n");
    break;
  default:
    return EXCEPTION_CONTINUE_SEARCH;
  }
  _clearfp();
  return EXCEPTION_CONTINUE_EXECUTION;
}
#else
static void fpe_signal_handler(int sig, siginfo_t *sip, void *scp) {
  // Use only async-signal-safe functions (write(2), not std::cerr/fprintf)
  static constexpr char prefix[] = "FPE (continuing): ";
  static constexpr char msg_div[] = "division by zero\n";
  static constexpr char msg_inv[] = "invalid operation\n";
  static constexpr char msg_ovf[] = "overflow\n";
  static constexpr char msg_unk[] = "unknown\n";
  write(STDERR_FILENO, prefix, sizeof(prefix) - 1);
  switch (sip->si_code) {
  case FPE_FLTDIV:
    write(STDERR_FILENO, msg_div, sizeof(msg_div) - 1);
    break;
  case FPE_FLTINV:
    write(STDERR_FILENO, msg_inv, sizeof(msg_inv) - 1);
    break;
  case FPE_FLTOVF:
    write(STDERR_FILENO, msg_ovf, sizeof(msg_ovf) - 1);
    break;
  default:
    write(STDERR_FILENO, msg_unk, sizeof(msg_unk) - 1);
    break;
  }

  // Clear FP exception flags in the saved signal context so the kernel
  // does not re-raise on sigreturn, then clear hardware state as well.
#if defined(__linux__) && (defined(__x86_64__) || defined(__i386__))
  ucontext_t *ctx = static_cast<ucontext_t *>(scp);
  if (ctx->uc_mcontext.fpregs) {
    ctx->uc_mcontext.fpregs->swd &= ~0x3Fu;   // x87 status word
    ctx->uc_mcontext.fpregs->mxcsr &= ~0x3Fu; // SSE MXCSR
  }
#endif
  feclearexcept(FE_ALL_EXCEPT);
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
  fprintf(stderr, "FPE trapping not supported on this platform.\n");
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
