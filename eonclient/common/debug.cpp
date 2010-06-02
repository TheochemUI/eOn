#include <stdarg.h>
#include <stdio.h>
#include "debug.h"
#if 0
void dbg(char const * fmt, ...)
{
    if(debug == 1)
    {
        va_list args;
        va_start(args, fmt);
        fprintf(stdout, fmt, args);
        va_end(args);
    }
}
#endif
