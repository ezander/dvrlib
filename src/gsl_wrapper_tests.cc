#include <cmath>

#include "gsl_wrapper.h"
#include "dvr_assert.h"

using namespace dvrlib;

void test_vector() {
  vector v(3, 4.0);
  dvr_assert(v.get(0) == 4.0);
  dvr_assert(v.get(1) == 4.0);
  dvr_assert(v.get(2) == 4.0);

  for(int i = 0; i < 3; i++) {
    v.set(i, 1.23 + i);
  }
  dvr_assert(v.get(0) == 1.23);
  dvr_assert(v.get(1) == 2.23);
  dvr_assert(v.get(2) == 3.23);

  vector w(3, 5.0);
  dvr_assert(w.get(0) == 5);

  w += v;
  dvr_assert(w.get(0) == 6.23);
  dvr_assert(w.get(2) == 8.23);
  dvr_assert(v.get(0) == 1.23);

  v = w * 3;
  dvr_assert(v.get(0) == 18.69);
  dvr_assert(w.get(0) == 6.23);

  double src[] = {3, 4, 10, 7};
  vector x(4, src);
  dvr_assert(x.get(0) == 3);
  dvr_assert(x.get(3) == 7);

  x -= x;
  dvr_assert(x.get(0) == 0);
  dvr_assert(x.get(1) == 0);

  vector z(4, src);
  dvr_assert(z.get(0) == 3);
  dvr_assert(z.get(3) == 7);

  x = -z;
  dvr_assert(x.get(0) == -3);
  dvr_assert(x.get(3) == -7);

  dvr_assert(x[0] == -3);
  dvr_assert(x[3] == -7);
}

void test_vector_view() {
  vector v(10, 4.0);
  vector_view vv = v.subvector(3, 2);
  vv.set(0, 2);
  vv.set(1, 5);
  dvr_assert(v.get(0) == 4);
  dvr_assert(v.get(3) == 2);
  dvr_assert(v.get(4) == 5);
  vector x(vv);
  x.set(0, 7);
  dvr_assert(x.get(0) == 7);
  dvr_assert(x.get(1) == 5);
  dvr_assert(vv.get(0) == 2);
  dvr_assert(v.get(3) == 2);
  x.set(1, 8);
  vv = x;
  dvr_assert(v.get(3) == 7);
  dvr_assert(v.get(4) == 8);
  x = vv;
  vv.set(0, -10);
  dvr_assert(x.get(0) == 7);
  dvr_assert(vv.get(0) == -10);

  vector v3(3);
  v3.set(0, 3);
  v3.set(1, -4);
  v3.set(2, 12);
  dvr_assert(v3.norm1() == 19);
  dvr_assert(v3.norm2() == 13);
}

void test_matrix() {
  matrix m(5, 6, true);
  dvr_assert(m.get(0, 0) == 1);
  dvr_assert(m.get(0, 1) == 0);
  dvr_assert(m.get(4, 4) == 1);
  dvr_assert(m.get(4, 5) == 0);

  m.set(3, 4, 7);
  m.set(1, 3, 9);
  matrix m2 = m.transpose();
  dvr_assert(m.get(3, 4) == 7);
  dvr_assert(m.get(1, 3) == 9);
  dvr_assert(m2.get(4, 3) == 7);
  dvr_assert(m2.get(3, 1) == 9);

  matrix A(2, 2);
  A.set(0, 0, 4);
  A.set(0, 1, 6);
  A.set(1, 0, 10);
  A.set(1, 1, 3);
  vector b(2);
  b.set(0, 38);
  b.set(1, 59);
  vector c = A.linsolve(b);
  dvr_assert(c.get(0) == 5);
  dvr_assert(c.get(1) == 3);
  dvr_assert(b.get(1) == 59);
  dvr_assert(A.get(0, 0) == 4);
  dvr_assert(A.get(1, 1) == 3);

  matrix B(A.inverse());
  double det = -48;
  dvr_assert(fabs(B.get(0, 0) - +3.0 / det) < 1e-9);
  dvr_assert(fabs(B.get(0, 1) - -6.0 / det) < 1e-9);
  dvr_assert(fabs(B.get(1, 0) - -10.0 / det) < 1e-9);
  dvr_assert(fabs(B.get(1, 1) - +4.0 / det) < 1e-9);

  matrix A2(2, 3);
  A2.set(0, 0, 4);
  A2.set(0, 1, 6);
  A2.set(0, 2, 7);
  A2.set(1, 0, 10);
  A2.set(1, 1, 3);
  A2.set(1, 2, 5);
  vector x(3);
  x.set(0, 1);
  x.set(1, 2);
  x.set(2, 3);
  vector y = A2 * x;
  dvr_assert(y.get(0) == 37);
  dvr_assert(y.get(1) == 31);

  double src[][3] = {{3, 4, 10}, {7, 5, 11}};
  matrix A3(2, 3, src);
  dvr_assert(A3.get(0, 0) == 3);
  dvr_assert(A3.get(1, 0) == 7);
  dvr_assert(A3.get(0, 1) == 4);
  dvr_assert(A3.get(1, 2) == 11);

  matrix A4(2, 3, src);
  A4 = A3 + A3;
  dvr_assert(A4.get(0, 0) == 6);
  dvr_assert(A4.get(1, 0) == 14);
  dvr_assert(A4.get(0, 1) == 8);
  dvr_assert(A4.get(1, 2) == 22);

  A4 += A3;
  dvr_assert(A4.get(0, 0) == 9);
  dvr_assert(A4.get(1, 0) == 21);
  dvr_assert(A4.get(0, 1) == 12);
  dvr_assert(A4.get(1, 2) == 33);

  A4 -= A3;
  dvr_assert(A4.get(0, 0) == 6);
  dvr_assert(A4.get(1, 0) == 14);
  dvr_assert(A4.get(0, 1) == 8);
  dvr_assert(A4.get(1, 2) == 22);

  A4 = A4 - A3;
  dvr_assert(A3.get(0, 0) == 3);
  dvr_assert(A3.get(1, 0) == 7);
  dvr_assert(A3.get(0, 1) == 4);
  dvr_assert(A3.get(1, 2) == 11);

  A4 = -A3;
  dvr_assert(A3.get(0, 0) == 3);
  dvr_assert(A3.get(1, 0) == 7);
  dvr_assert(A3.get(0, 1) == 4);
  dvr_assert(A3.get(1, 2) == 11);

  dvr_assert(A4.get(0, 0) == -3);
  dvr_assert(A4.get(1, 0) == -7);
  dvr_assert(A4.get(0, 1) == -4);
  dvr_assert(A4.get(1, 2) == -11);

  A4 *= -2.0;
  dvr_assert(A4.get(0, 0) == 6);
  dvr_assert(A4.get(1, 0) == 14);
  dvr_assert(A4.get(0, 1) == 8);
  dvr_assert(A4.get(1, 2) == 22);

  A4 = A4 * -2.0;
  dvr_assert(A4.get(0, 0) == -12);
  dvr_assert(A4.get(1, 0) == -28);
  dvr_assert(A4.get(0, 1) == -16);
  dvr_assert(A4.get(1, 2) == -44);

  vector d(A4[1]);
  dvr_assert(d.get(0) == -28);
  dvr_assert(d.get(1) == -20);
  dvr_assert(d.get(2) == -44);
}

void test_matrix_view() {
  matrix m(10, 10, true);
  matrix_view mv = m.submatrix(3, 3, 2, 2);
  dvr_assert(mv.get(0, 0) == 1);
  dvr_assert(mv.get(0, 1) == 0);

  mv.set(0, 1, 2);
  mv.set(1, 0, 5);
  dvr_assert(mv.get(0, 1) == 2);
  dvr_assert(mv.get(1, 0) == 5);
  dvr_assert(m.get(3, 4) == 2);
  dvr_assert(m.get(4, 3) == 5);

  matrix m3(2, 2, true);
  m3.set(0, 0, 5);
  m3.set(1, 0, 6);
  m3.set(0, 1, 7);
  m3.set(1, 1, 8);
  m.submatrix(2, 4, 2, 2) = m3;
  dvr_assert(m.get(2, 4) == 5);
  dvr_assert(m.get(3, 4) == 6);
  dvr_assert(m.get(2, 5) == 7);
  dvr_assert(m.get(3, 5) == 8);
}

void gsl_wrapper_test_suite() {
  test_vector();
  test_vector_view();
  test_matrix();
  test_matrix_view();
}
