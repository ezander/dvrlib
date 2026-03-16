#include "gsl_wrapper.h"
#include "vdi2048.h"

using namespace dvrlib;

int main() {
  gsl_enable_exceptions();
  example_VDI2048();
  return 0;
}
