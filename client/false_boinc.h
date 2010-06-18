#ifndef FALSE_BOINC_H
#define FALSE_BOINC_H
#include <cstdlib>
#include <cstring>
// For degugging
#define _BOINC_API_
#define _BOINC_DIAGNOSTICS_
#define _FILESYS_

inline int boinc_init() {return 0;}
inline void boinc_finish(int) {}
inline int boinc_resolve_filename(const char *str1, char *str2, size_t n)
{
      strncpy(str2, str1, n);
      return 0;
}

#define boinc_fopen fopen

#endif
