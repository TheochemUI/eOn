//Includes for FPE trapping

/* GH: fix for arm
#ifdef OSX
    #include <xmmintrin.h>
    #include <mach/mach_init.h>
    #include <mach/task.h>
#endif
*/
#if defined(__APPLE__) && defined(__x86_64__)
    #include <xmmintrin.h>
#endif

#ifdef LINUX
    #include <fenv.h>
#endif

#include "ExceptionsEON.h"

#include <cfenv>

void enableFPE(void)
{
    // Floating Point Trapping. It is platform specific!
    // This causes the program to crash on divison by zero,
    // invalid operations, and overflows.

/*
    #ifdef LINUX
        feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);
    #endif

    #ifdef OSX
        _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK()
                               & ~_MM_MASK_INVALID 
                               & ~_MM_MASK_DIV_ZERO
                               & ~_MM_MASK_OVERFLOW);
    #endif 
*/

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
}

void disableFPE(void)
{
    #ifdef LINUX
        fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    #endif

// GH: this all needs to be fixed for x86_64 vs aarch64
/*
    #ifdef OSX
        _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK()
                               | _MM_MASK_INVALID 
                               | _MM_MASK_DIV_ZERO
                               | _MM_MASK_OVERFLOW);
    #endif
*/

    #if defined(__APPLE__) && defined(__x86_64__)
        // Enable floating-point exceptions on Intel macOS
        _MM_SET_EXCEPTION_MASK(
            _MM_MASK_MASK &
            ~(_MM_MASK_INVALID | _MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW));
    #endif
}

int getFPEState(void) {
    #ifdef LINUX
        return fegetexcept();
    #endif

// GH: this needs to be fixed for x86_64 vs aarch64
/*
    #ifdef OSX
        return _MM_GET_EXCEPTION_MASK();
    #endif
*/
    #if defined(__APPLE__) && defined(__x86_64__)
        return _MM_GET_EXCEPTION_MASK();
    #endif
}

bool isFPEEnabled(void) {
    int currentFPEState = getFPEState();
    #ifdef LINUX
        return (currentFPEState & (FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW)) != 0;
    #endif
// GH: this needs to be fixed for x86_64 vs aarch64
/*
    #ifdef OSX
        return (currentFPEState & (_MM_MASK_INVALID | _MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW)) != 0;
    #endif
*/
    #if defined(__APPLE__) && defined(__x86_64__)
        return (currentFPEState & (_MM_MASK_INVALID | _MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW)) != 0;
    #endif

    return false;  // Fallback, should never reach here
}
