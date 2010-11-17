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
#include <csignal>
#include <iostream>
#define debug 1
//void dbg(char const *, ...);
#define dbg(x) error
#define _file printf("%s:%d\n", __FILE__, __LINE__);
#define _fun printf("%s:%d\n", __FUNCTION__, __LINE__);
#define printl(x) std::cout << # x << " = " << x << std::endl;
#if 0
      #undef assert
      #define assert(x) \
            if (not (x)) { \
                  std::cerr << # x << std::endl; \
                  std::cerr << __FILE__ << ':' << __FUNCTION__ << ':' << __LINE__ << std::endl; \
                   raise(SIGINT); \
            }
#endif
