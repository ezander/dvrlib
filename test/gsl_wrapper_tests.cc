#include <catch2/catch_all.hpp>
#include <cmath>
#include "gsl_wrapper.h"

using namespace dvrlib;

TEST_CASE("vector") {
    vector v(3, 4.0);
    REQUIRE(v.get(0) == 4.0);
    REQUIRE(v.get(1) == 4.0);
    REQUIRE(v.get(2) == 4.0);

    for (int i = 0; i < 3; i++) {
        v.set(i, 1.23 + i);
    }
    REQUIRE(v.get(0) == 1.23);
    REQUIRE(v.get(1) == 2.23);
    REQUIRE(v.get(2) == 3.23);

    vector w(3, 5.0);
    REQUIRE(w.get(0) == 5);

    w += v;
    REQUIRE(w.get(0) == 6.23);
    REQUIRE(w.get(2) == 8.23);
    REQUIRE(v.get(0) == 1.23);

    v = w * 3;
    REQUIRE(v.get(0) == 18.69);
    REQUIRE(w.get(0) == 6.23);

    double src[] = {3, 4, 10, 7};
    vector x(4, src);
    REQUIRE(x.get(0) == 3);
    REQUIRE(x.get(3) == 7);

    x -= x;
    REQUIRE(x.get(0) == 0);
    REQUIRE(x.get(1) == 0);

    vector z(4, src);
    REQUIRE(z.get(0) == 3);
    REQUIRE(z.get(3) == 7);

    x = -z;
    REQUIRE(x.get(0) == -3);
    REQUIRE(x.get(3) == -7);

    REQUIRE(x[0] == -3);
    REQUIRE(x[3] == -7);
}

TEST_CASE("vector_view") {
    vector v(10, 4.0);
    vector_view vv = v.subvector(3, 2);
    vv.set(0, 2);
    vv.set(1, 5);
    REQUIRE(v.get(0) == 4);
    REQUIRE(v.get(3) == 2);
    REQUIRE(v.get(4) == 5);

    vector x(vv);
    x.set(0, 7);
    REQUIRE(x.get(0) == 7);
    REQUIRE(x.get(1) == 5);
    REQUIRE(vv.get(0) == 2);
    REQUIRE(v.get(3) == 2);

    x.set(1, 8);
    vv = x;
    REQUIRE(v.get(3) == 7);
    REQUIRE(v.get(4) == 8);

    x = vv;
    vv.set(0, -10);
    REQUIRE(x.get(0) == 7);
    REQUIRE(vv.get(0) == -10);

    vector v3(3);
    v3.set(0, 3);
    v3.set(1, -4);
    v3.set(2, 12);
    REQUIRE(v3.norm1() == 19);
    REQUIRE(v3.norm2() == 13);
}

TEST_CASE("matrix") {
    matrix m(5, 6, true);
    REQUIRE(m.get(0, 0) == 1);
    REQUIRE(m.get(0, 1) == 0);
    REQUIRE(m.get(4, 4) == 1);
    REQUIRE(m.get(4, 5) == 0);

    m.set(3, 4, 7);
    m.set(1, 3, 9);
    matrix m2 = m.transpose();
    REQUIRE(m.get(3, 4) == 7);
    REQUIRE(m.get(1, 3) == 9);
    REQUIRE(m2.get(4, 3) == 7);
    REQUIRE(m2.get(3, 1) == 9);

    matrix A(2, 2);
    A.set(0, 0, 4);
    A.set(0, 1, 6);
    A.set(1, 0, 10);
    A.set(1, 1, 3);
    vector b(2);
    b.set(0, 38);
    b.set(1, 59);
    vector c = A.linsolve(b);
    REQUIRE(c.get(0) == 5);
    REQUIRE(c.get(1) == 3);
    REQUIRE(b.get(1) == 59);
    REQUIRE(A.get(0, 0) == 4);
    REQUIRE(A.get(1, 1) == 3);

    matrix B(A.inverse());
    double det = -48;
    REQUIRE(fabs(B.get(0, 0) - +3.0 / det) < 1e-9);
    REQUIRE(fabs(B.get(0, 1) - -6.0 / det) < 1e-9);
    REQUIRE(fabs(B.get(1, 0) - -10.0 / det) < 1e-9);
    REQUIRE(fabs(B.get(1, 1) - +4.0 / det) < 1e-9);

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
    REQUIRE(y.get(0) == 37);
    REQUIRE(y.get(1) == 31);

    double src[][3] = {{3, 4, 10}, {7, 5, 11}};
    matrix A3(2, 3, src);
    REQUIRE(A3.get(0, 0) == 3);
    REQUIRE(A3.get(1, 0) == 7);
    REQUIRE(A3.get(0, 1) == 4);
    REQUIRE(A3.get(1, 2) == 11);

    matrix A4(2, 3, src);
    A4 = A3 + A3;
    REQUIRE(A4.get(0, 0) == 6);
    REQUIRE(A4.get(1, 0) == 14);
    REQUIRE(A4.get(0, 1) == 8);
    REQUIRE(A4.get(1, 2) == 22);

    A4 += A3;
    REQUIRE(A4.get(0, 0) == 9);
    REQUIRE(A4.get(1, 0) == 21);
    REQUIRE(A4.get(0, 1) == 12);
    REQUIRE(A4.get(1, 2) == 33);

    A4 -= A3;
    REQUIRE(A4.get(0, 0) == 6);
    REQUIRE(A4.get(1, 0) == 14);
    REQUIRE(A4.get(0, 1) == 8);
    REQUIRE(A4.get(1, 2) == 22);

    A4 = A4 - A3;
    REQUIRE(A3.get(0, 0) == 3);
    REQUIRE(A3.get(1, 0) == 7);
    REQUIRE(A3.get(0, 1) == 4);
    REQUIRE(A3.get(1, 2) == 11);

    A4 = -A3;
    REQUIRE(A3.get(0, 0) == 3);
    REQUIRE(A3.get(1, 0) == 7);
    REQUIRE(A3.get(0, 1) == 4);
    REQUIRE(A3.get(1, 2) == 11);
    REQUIRE(A4.get(0, 0) == -3);
    REQUIRE(A4.get(1, 0) == -7);
    REQUIRE(A4.get(0, 1) == -4);
    REQUIRE(A4.get(1, 2) == -11);

    A4 *= -2.0;
    REQUIRE(A4.get(0, 0) == 6);
    REQUIRE(A4.get(1, 0) == 14);
    REQUIRE(A4.get(0, 1) == 8);
    REQUIRE(A4.get(1, 2) == 22);

    A4 = A4 * -2.0;
    REQUIRE(A4.get(0, 0) == -12);
    REQUIRE(A4.get(1, 0) == -28);
    REQUIRE(A4.get(0, 1) == -16);
    REQUIRE(A4.get(1, 2) == -44);

    vector d(A4[1]);
    REQUIRE(d.get(0) == -28);
    REQUIRE(d.get(1) == -20);
    REQUIRE(d.get(2) == -44);
}

TEST_CASE("matrix_view") {
    matrix m(10, 10, true);
    matrix_view mv = m.submatrix(3, 3, 2, 2);
    REQUIRE(mv.get(0, 0) == 1);
    REQUIRE(mv.get(0, 1) == 0);

    mv.set(0, 1, 2);
    mv.set(1, 0, 5);
    REQUIRE(mv.get(0, 1) == 2);
    REQUIRE(mv.get(1, 0) == 5);
    REQUIRE(m.get(3, 4) == 2);
    REQUIRE(m.get(4, 3) == 5);

    matrix m3(2, 2, true);
    m3.set(0, 0, 5);
    m3.set(1, 0, 6);
    m3.set(0, 1, 7);
    m3.set(1, 1, 8);
    m.submatrix(2, 4, 2, 2) = m3;
    REQUIRE(m.get(2, 4) == 5);
    REQUIRE(m.get(3, 4) == 6);
    REQUIRE(m.get(2, 5) == 7);
    REQUIRE(m.get(3, 5) == 8);
}
