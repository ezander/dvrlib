#ifndef DVRLIB_DVR_ASSERT_H
#define DVRLIB_DVR_ASSERT_H

#include <stdexcept>
#include <string>

namespace dvrlib {

struct assertion_error : public std::logic_error {
  assertion_error(const char* cond, const char* file, int line)
      : std::logic_error(std::string("assertion failed: ") + cond +
                         " (" + file + ":" + std::to_string(line) + ")") {}
};

}  // namespace dvrlib

#define dvr_assert(cond) \
  do { \
    if (!(cond)) \
      throw dvrlib::assertion_error(#cond, __FILE__, __LINE__); \
  } while (0)

#endif  // DVRLIB_DVR_ASSERT_H
