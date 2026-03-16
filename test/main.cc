#include <catch2/catch_all.hpp>
#include "gsl_wrapper.h"

using namespace dvrlib;

int main(int argc, char* argv[]) {
    gsl_enable_exceptions();
    return Catch::Session().run(argc, argv);
}
