cmake_minimum_required(VERSION 3.15)
project(dvrlib VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GSL REQUIRED)

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

# Test / example executable
add_executable(dvrlib_main
    src/main.cc
    src/gsl_wrapper_tests.cc
    src/recon_tests.cc
)
target_link_libraries(dvrlib_main PRIVATE dvrlib)
target_compile_options(dvrlib_main PRIVATE -Wall)

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
