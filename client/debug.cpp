//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//
//-----------------------------------------------------------------------------------
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
