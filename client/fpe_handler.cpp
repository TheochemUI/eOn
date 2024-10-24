#include "fpe_handler.h"

#include <cfenv>
#include <cmath>
#include <csignal>
#include <iostream>

#ifdef _WIN32
#include <float.h>
#endif

#if defined(__APPLE__) && defined(__x86_64__)
#include <xmmintrin.h>
#endif

namespace eonc {

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

void enableFPE() {
#ifdef _WIN32
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

  // Register signal handler
  struct sigaction act;
  act.sa_sigaction = fpe_signal_handler;
  sigemptyset(&act.sa_mask);
  act.sa_flags = SA_SIGINFO;
  sigaction(SIGFPE, &act, nullptr);
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
