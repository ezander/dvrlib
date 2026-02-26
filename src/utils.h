#ifndef DVRLIB_UTILS_H
#define DVRLIB_UTILS_H

#define PRINT(x) std::cout << #x << ": " << (x) << std::endl
#define PRINT_SEP std::cout << "======================================================================" << std::endl
#define PRINT_TITLE(x)                      \
  PRINT_SEP;                                \
  std::cout << "     " << (x) << std::endl; \
  PRINT_SEP

#endif  // DVRLIB_UTILS_H
