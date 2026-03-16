#include <iostream>
#include <vector>

#include "gsl_wrapper.h"
#include "recon.h"
#include "utils.h"
#include "recon_system.h"
#include "vdi2048.h"

using namespace dvrlib;

void example_VDI2048_1(bool increase_var = false, bool keep_free = true) {
  recon_system system;

  system.add_var("m_FDKeI", 46.241, 0.800);
  system.add_var("m_FDKeII", 45.668, 0.790);
  system.add_var("m_SpI", 44.575, 0.535);
  system.add_var("m_SpII", 44.319, 0.532);
  system.add_var("m_V", 0.525, 0.105);
  system.add_var("m_HK", 69.978, 0.854);
  system.add_var("m_A7", 10.364, 0.168);
  system.add_var("m_A6", 3.744, 0.058);
  system.add_var("m_A5", 4.391, 0.058);
  system.add_var("m_HDNK", 18.498, 0.205);
  system.add_var("m_D", 2.092, 0.272);

  if(keep_free) {
    system.add_var("m_FD1", 0, -1);
    system.add_var("m_FD2", 0, -1);
    system.add_var("m_FD3", 0, -1);
    system.add_var("m_HDAnz", 0, -1);
  }

  system.add_covariance_coeff("m_FDKeI", "m_FDKeII", 0.2);
  system.add_covariance_coeff("m_SpI", "m_SpII", 0.4);

  if(increase_var) {
    system.change_var("m_FDKeI", 46.241, 2.500);
    system.change_var("m_FDKeII", 45.668, 2.500);
  }

  // get and print the covariance matrix
  matrix S_x = system.get_covariance_matrix();
  PRINT_TITLE("Matrix S_X (compare A9, values from A6 and A8)");
  PRINT(S_x);

  // define the constraints
  //   m_FD1 = m_FDKeI + m_FDKeII - 0.2 * m_V
  //   m_FD2 = m_SpI + m_SpII - 0.6 * m_V
  //   m_FD3 = m_HK + m_A7 + m_A6 + m_A5 + 0.4 * m_V
  //   m_HDAnz = m_A7 + m_A6 + m_A5
  // plus (if keep_free == true)
  //   M_FD1 - M_FD2 = 0
  //   M_FD2 - M_FD3 = 0
  //   M_HDAnz - M_HDNK = 0
  auto F_with_free = []() -> matrix {
    double Fc[][15] = {
      {1, 1, 0, 0, -0.2, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0},
      {0, 0, 1, 1, -0.6, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0},
      {0, 0, 0, 0, 0.4, 1, 1, 1, 1, 0, 0, 0, 0, -1, 0},
      {0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, -1},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1},
    };
    return matrix(7, 15, Fc);
  };

  auto F_without_free = []() -> matrix {
    double Fc[][11] = {
      {1, 1, -1, -1, 0.4, 0, 0, 0, 0, 0, 0},
      {0, 0, 1, 1, -1, -1, -1, -1, -1, 0, 0},
      {0, 0, 0, 0, 0, 0, 1, 1, 1, -1, 0},
    };
    return matrix(3, 11, Fc);
  };

  matrix F = keep_free ? F_with_free() : F_without_free();

  PRINT_TITLE("Constraints (compare A1-A4 and A12)");
  system.print_constraints(F);

  PRINT_TITLE("Constraints as Matrix F + (+extra cols)");
  PRINT(F);

  vector x = system.get_values();

  // A15
  // PRINT_TITLE("Matrix F*S_X (A15)");
  // PRINT(F*S_x);

  PRINT_TITLE("Matrix F*X (A17)");
  PRINT(F * x);

  PRINT_SEP;

  vector v(x.size());

  lin_recon(F * x, S_x, F, v);

  PRINT_TITLE("Vector v (compare A20)");
  PRINT(v);

  // test that the constraints are fulfilled
  PRINT_TITLE("Constraint fulfillment F*(v+x) (should be close to zero)");
  PRINT(F * (x + v));

  // compute and print covariance S_v
  matrix S_v(S_x.size1(), S_x.size2());
  lin_cov_update(S_x, F, S_v);

  PRINT_TITLE("Matrix S_v (compare A21)");
  PRINT(S_v);

  // compute and print covariance update S_X_new
  matrix S_x_new(S_x.size1(), S_x.size2());
  S_x_new = S_x - S_v;

  PRINT_TITLE("Matrix S_x_new (compare A33)");
  PRINT(S_x_new);

  vector x_new(x.size());
  x_new = x + v;

  PRINT_TITLE("Vector x_new (compare A34)");
  PRINT(x_new);

  vector confint_new(S_x.size1());
  extract_confidence(S_x_new, confint_new);

  PRINT_TITLE("confint_new (compare A34)");
  PRINT(confint_new);

  PRINT_TITLE("Original variables and confidences");
  system.print_vars();

  PRINT_TITLE("Updated variables and confidences");
  auto updated_sys = system.updated(x_new, confint_new);
  updated_sys.print_vars();

}

void example_VDI2048() {
  example_VDI2048_1(true, true);
}