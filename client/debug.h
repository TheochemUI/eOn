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
