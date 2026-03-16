#ifndef DVRLIB_UTILS_H
#define DVRLIB_UTILS_H

#include <cstdio>
#include <sstream>
#include <string>

#define PRINT(x)                           \
  do {                                     \
    std::ostringstream _s;                 \
    _s << (x);                             \
    printf(#x ": %s\n", _s.str().c_str()); \
  } while(0)
#define PRINT_SEP \
  printf("======================================================================\n")
#define PRINT_TITLE(x)                        \
  do {                                        \
    PRINT_SEP;                                \
    printf("     %s\n", (std::string() + (x)).c_str()); \
    PRINT_SEP;                                \
  } while(0)

#endif  // DVRLIB_UTILS_H
