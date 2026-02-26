# dvrlib

A C++ library for **Data Validation and Reconciliation** (DVR) based on the methods described in the [VDI 2048](https://www.vdi.de/richtlinien/details/vdi-2048-blatt-1-messunsicherheiten-bei-abnahmemessungen-an-energie-und-kraftwerkstechnischen-anlagen-grundlagen) standard.

Given a set of measured values with known uncertainties and a system of constraints (e.g. mass or energy balances), dvrlib computes reconciled values that satisfy the constraints while minimising the weighted corrections according to the measurement covariances. It also propagates uncertainties to produce updated confidence intervals for the reconciled results.

## Features

- RAII C++ wrappers around GSL vectors and matrices with operator overloading
- Linear reconciliation (`lin_recon`) and covariance update (`lin_cov_update`)
- Iterative non-linear reconciliation via linearisation (`recon`)
- High-level `recon_system` class for defining variables with names, measured values, confidence intervals, and covariance coefficients
- Conversion between 95% confidence intervals and variances (`confint2var`, `var2confint`)
- Worked example reproducing VDI 2048 Appendix A

## Prerequisites

- A C++ compiler with C++11 support (e.g. GCC)
- [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
- CMake 3.15+

On Debian/Ubuntu:

```bash
sudo apt install libgsl-dev cmake build-essential
```

## Building

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug   # configure
cmake --build build                        # build library + test binary
cmake --build build --target doc           # generate API docs (requires Doxygen)
rm -rf build                               # clean
```

Produces `build/libdvrlib.a` (static library) and `build/dvrlib_main` (test binary).

## Running Tests

```bash
cmake --build build && build/dvrlib_main
```

The test binary runs the GSL wrapper tests, the reconciliation tests, and the VDI 2048 example. A clean exit (return code 0) means all tests pass; an assertion failure or uncaught exception indicates a failure.

## Usage Example

```cpp
#include "gsl_wrapper.h"
#include "recon.h"
#include "recon_system.h"

using namespace dvrlib;

int main() {
    gsl_enable_exceptions();

    // Define the measurement system
    recon_system system;
    system.add_var("flow_in",  100.0, 2.0);   // measured value, 95% confidence interval
    system.add_var("flow_out",  98.5, 1.5);
    system.add_var("loss",       0.8, 0.5);

    // Constraint: flow_in - flow_out - loss = 0
    double Fc[][3] = {{1, -1, -1}};
    matrix F(1, 3, Fc);

    // Reconcile
    matrix S_x = system.get_covariance_matrix();
    vector x = system.get_values();
    vector v(x.size());
    lin_recon(F * x, S_x, F, v);

    // x + v now satisfies the constraint within measurement uncertainties
}
```

## Project Structure

```
src/
  gsl_wrapper.{h,cc}       GSL vector/matrix RAII wrappers
  gsl_wrapper_tests.{h,cc} Tests for the wrapper layer
  recon.{h,cc}              Reconciliation algorithms
  recon_system.{h,cc}       High-level system description
  recon_tests.{h,cc}        Tests for reconciliation
  vdi2048.{h,cc}            VDI 2048 worked example
  main.cc                   Test runner entry point
```

## License

MIT License. See [LICENSE](LICENSE) for details.
