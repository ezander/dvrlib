cmake_minimum_required(VERSION 3.15)
project(dvrlib VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GSL REQUIRED)

include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.8.0
)
FetchContent_MakeAvailable(Catch2)

# Core library
add_library(dvrlib STATIC
    src/gsl_wrapper.cc
    src/recon.cc
    src/recon_system.cc
    src/vdi2048.cc
)
target_include_directories(dvrlib PUBLIC src)
target_link_libraries(dvrlib PUBLIC GSL::gsl GSL::gslcblas)
target_compile_options(dvrlib PRIVATE -Wall)

# Tests
add_executable(dvrlib_main
    test/main.cc
    test/gsl_wrapper_tests.cc
    test/recon_tests.cc
    test/vdi2048_test.cc
)
target_link_libraries(dvrlib_main PRIVATE dvrlib Catch2::Catch2)
target_compile_options(dvrlib_main PRIVATE -Wall)

include(CTest)
include(Catch)
catch_discover_tests(dvrlib_main)

# Demo
add_executable(dvrlib_demo demo/vdi2048_demo.cc)
target_link_libraries(dvrlib_demo PRIVATE dvrlib)
target_compile_options(dvrlib_demo PRIVATE -Wall)

# Documentation (optional — requires Doxygen)
find_package(Doxygen)
if(DOXYGEN_FOUND)
    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating documentation with Doxygen"
        VERBATIM
    )
endif()
