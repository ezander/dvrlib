#include <catch2/catch_all.hpp>
#include "gsl_wrapper.h"
#include "vdi2048.h"

using namespace dvrlib;

int main(int argc, char* argv[]) {
    gsl_enable_exceptions();
    return Catch::Session().run(argc, argv);
}

TEST_CASE("VDI2048 example") {
    example_VDI2048();
}
