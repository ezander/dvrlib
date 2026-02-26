#include <cstdio>
#include <gsl/gsl_errno.h>
#include "gsl_wrapper.h"
#include "gsl_wrapper_tests.h"
#include "recon_tests.h"
#include "vdi2048.h"

using namespace dvrlib;

int main(void) {
  gsl_enable_exceptions();

  try {
    gsl_wrapper_test_suite();
    recon_test_suite();
    example_VDI2048();
  } catch(const gsl_exception& e) {
    printf("Caught GSL exception\n");
    printf("  reason:    %s\n", e.reason);
    printf("  file:      %s\n", e.file);
    printf("  line:      %d\n", e.line);
    printf("  gsl_errno: %d\n", e.gsl_errno);
  }

  return 0;
}
