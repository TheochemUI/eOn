//Includes for FPE trapping
#ifdef OSX
    #include <xmmintrin.h>
    #include <mach/mach_init.h>
    #include <mach/task.h>
#endif
#ifdef LINUX
    #include <fenv.h>
#endif

#include "ExceptionsEON.h"


void enableFPE(void)
{
    // Floating Point Trapping. It is platform specific!
    // This causes the program to crash on divison by zero,
    // invalid operations, and overflows.
    #ifdef LINUX
        feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);
    #endif
    #ifdef OSX
        _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK()
                               & ~_MM_MASK_INVALID 
                               & ~_MM_MASK_DIV_ZERO
                               & ~_MM_MASK_OVERFLOW);
    #endif 
}

void disableFPE(void)
{
    #ifdef LINUX
        fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    #endif
    #ifdef OSX
        _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK()
                               | _MM_MASK_INVALID 
                               | _MM_MASK_DIV_ZERO
                               | _MM_MASK_OVERFLOW);
    #endif 
}

int getFPEState(void) {
    #ifdef LINUX
        return fegetexcept();
    #endif
    #ifdef OSX
        return _MM_GET_EXCEPTION_MASK();
    #endif
}

bool isFPEEnabled(void) {
    int currentFPEState = getFPEState();
    #ifdef LINUX
        return (currentFPEState & (FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW)) != 0;
    #endif
    #ifdef OSX
        return (currentFPEState & (_MM_MASK_INVALID | _MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW)) != 0;
    #endif
    return false;  // Fallback, should never reach here
}