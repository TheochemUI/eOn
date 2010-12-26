//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef FALSEBOINC_H
#define FALSEBOINC_H

#include <cstring>
#include <stdlib.h>
// For degugging
#define _BOINC_API_
#define _BOINC_DIAGNOSTICS_
#define _FILESYS_

inline int boinc_init() {return 0;}
inline void boinc_finish(int status) {exit(status);}
inline int boinc_resolve_filename(const char *str1, char *str2, size_t n)
{
      strncpy(str2, str1, n);
      return 0;
}

inline void boinc_fraction_done(double frac)
{
    #ifdef BOINC_DEBUG
    printf("%.1f%% done\n", frac*100);
    #endif
}
#define boinc_fopen fopen

#endif
