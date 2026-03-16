#include "gsl_wrapper.h"
#include "vdi2048.h"

#include <cstdio>
#include <cstring>

using namespace dvrlib;

int main(int argc, char* argv[]) {
  bool keep_free = false;

  for(int i = 1; i < argc; i++) {
    if(strcmp(argv[i], "--keep-free") == 0) {
      keep_free = true;
    } else if(strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
      printf("Usage: %s [options]\n", argv[0]);
      printf("\n");
      printf("Options:\n");
      printf("  --keep-free   Keep free variables in the system of equations\n");
      printf("                (Z-algorithm). Default: eliminate free variables first.\n");
      printf("  -h, --help    Show this help message.\n");
      return 0;
    } else {
      printf("Unknown option: %s\n", argv[i]);
      printf("Run with --help for usage.\n");
      return 1;
    }
  }

  gsl_enable_exceptions();
  example_VDI2048(keep_free);
  return 0;
}
