#include "fpe_handler.h"

#include <cfenv>
#include <iostream>

#ifdef _WIN32
#include <float.h>
#endif

void enableFPE() {
#ifdef _WIN32
  // Enable floating-point exceptions on Windows
  _controlfp_s(nullptr, 0, _MCW_EM);
  _controlfp_s(nullptr, ~(_EM_ZERODIVIDE | _EM_INVALID | _EM_OVERFLOW),
               _MCW_EM);
#elif defined(__unix__) || (defined(__APPLE__) && defined(__aarch64__))
  // Enable floating-point exceptions on Unixes and M1 MacOS
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#elif defined(__x86_64__)
  // Enable floating-point exceptions on Intel macOS
  fenv_t fenv;
  fegetenv(&fenv);
  fenv.__mxcsr &= ~(_MM_MASK_INVALID | _MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW);
  fesetenv(&fenv);
#else
  std::cerr << "FPE trapping not supported on this platform." << std::endl;
#endif
}
