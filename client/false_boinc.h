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
