#include <cassert>
#include <cmath>

#include "recon_tests.h"
#include "recon.h"
#include "gsl_wrapper.h"
#include "utils.h"

#define assert_equal(a, b)\
  assert((a)==(b))
#define assert_almost_equal(a, b)\
  assert(fabs((a)-(b))<1e-7*(fabs(a)+fabs(b)))
#define assert_vector_almost_equal_upto(a, b,eps)			\
  assert(((a)+(-1*(b))).norm2()<(eps)*((a).norm2()+(b).norm2()))
#define assert_vector_almost_equal(a, b)\
  assert(((a)+(-1*(b))).norm2()<1e-7*((a).norm2()+(b).norm2()))

namespace simple_system {
  // m1-m2=0
  // sigma1**2=1, sigma2**2=4
  double F[][2] = {{1, -1}};
  double S_x_diag[] = {1, 4};
  double S_x_inv_diag[] = {1, 0.25};
  double x[] = {0, 3};
  double v[] = {3.0*1.0/5.0, -3.0*4.0/5.0};
}

namespace simple_system2 {
  // m1-m3=0
  // m2-m4=0
  // m3-m4=0
  // sigma1**2=1, sigma3**2=4
  double F[][4] = {{1, 0, -1, 0}, {0, 1, 0, -1}, {0, 0, 1, -1}};
  double S_x_diag[] = {1, 4};
  double S_x_inv_diag[] = {1, 0.25};
  double x[] = {0, 3, 0, 0};
  double v[] = {3.0/5.0, -12.0/5.0, 3.0/5.0, 3.0/5.0};
}

namespace mischer_teiler {
  // m1+m2-m5=0
  // -m3-m4+m5=0
  double F[][5] = {{1, 1, 0, 0, -1}, {0, 0, -1, -1, 1}};
  double S_x_diag[] = {0.25, 0.25, 0.25, 0.25};
  double S_x_inv_diag[] = {4.0, 4.0, 4.0, 4.0};
  double x[] = {98, 99, 100, 100, 0};
  double v[] = {0.75, 0.75, -0.75, -0.75, 198.5};
}

void test_lin_recon() {
  {
    matrix F(2, 5, mischer_teiler::F);
    matrix S_x(4, 4, true, mischer_teiler::S_x_diag);
    vector x(5, mischer_teiler::x);

    vector r=F*x;
    vector v(5), v_exp(5, mischer_teiler::v);
  
    lin_recon(r, S_x, F, v);
    assert_vector_almost_equal(v, v_exp);
  }
  {
    matrix F(1, 2, simple_system::F);
    matrix S_x(2, 2, true, simple_system::S_x_diag);
    vector x(2, simple_system::x);

    vector r=F*x;
    vector v(2), v_exp(2, simple_system::v);
  
    lin_recon(r, S_x, F, v);
    assert_vector_almost_equal(v, v_exp);
  }
  {
    matrix F(3, 4, simple_system2::F);
    matrix S_x(2, 2, true, simple_system2::S_x_diag);
    vector x(4, simple_system2::x);

    vector r=F*x;
    vector v(4), v_exp(4, simple_system2::v);
  
    lin_recon(r, S_x, F, v);
    assert_vector_almost_equal(v, v_exp);
  }

}

void test_lin_recon_update() {
  matrix F(2, 5, mischer_teiler::F);
  matrix S_x_inv(4, 4, true, mischer_teiler::S_x_inv_diag);
  vector x(5, mischer_teiler::x);

  vector v(5), dv(5), v_exp(5, mischer_teiler::v);
  v.set(1, 3);
  v.set(2, 4);
  v.set(4, 7);
  vector r=F*(x+v);
  
  lin_recon_update(r, S_x_inv, F, v, dv);
  assert_vector_almost_equal(v+dv, v_exp);
}



class EnbiproDummy {
public:
  virtual vector getValues() = 0;
  virtual matrix getCovarianceMatrix() = 0;
  virtual matrix getJacobian(const vector& x) = 0;
  virtual vector getResidual(const vector& x) = 0;
};

class LinearEnbiproDummy : public EnbiproDummy {
public:
  vector getValues() {
    vector x(5, mischer_teiler::x);
    return x;
  }

  matrix getCovarianceMatrix() {
    matrix S_x(4, 4, true, mischer_teiler::S_x_diag);
    return S_x;
  }

  matrix getJacobian(const vector& x) {
    matrix F(2, 5, mischer_teiler::F);
    return F;
  }

  vector getResidual(const vector& x) {
    return getJacobian(x) * x;
  }
};


class QuadraticEnbiproDummy : public EnbiproDummy {
public:
  vector getValues() {
    vector x(5, mischer_teiler::x);
    x.set(4, 201);
    return x;
  }

  matrix getCovarianceMatrix() {
    matrix S_x(4, 4, true, mischer_teiler::S_x_diag);
    return S_x;
  }

  matrix getJacobian(const vector& x) {
    double d1 = 2*(x.get(0)+x.get(1)-x.get(4));
    double d2 = 2*(x.get(2)+x.get(3)-x.get(4));
    double J[][5] = {{d1, d1, 0, 0, -d1}, {0, 0, d2, d2, -d2}};
    matrix F(2, 5, J);
    return F;
  }

  vector getResidual(const vector& x) {
    vector r(2);
    r.set(0, pow( x.get(0)+x.get(1)-x.get(4),2));
    r.set(1, pow(-x.get(2)-x.get(3)+x.get(4),2));
    return r;
  }
};


class EnbiJacobianFunc : public func<vector,matrix> {
  EnbiproDummy* enbi;
public:
  EnbiJacobianFunc(EnbiproDummy* _enbi) :
    enbi(_enbi) {}
  virtual matrix operator()(const vector& arg);
};
matrix EnbiJacobianFunc::operator()(const vector& arg) {
  return enbi->getJacobian(arg);
}


class EnbiResidualFunc : public func<vector,vector> {
  EnbiproDummy* enbi;
public:
  EnbiResidualFunc(EnbiproDummy* _enbi) :
    enbi(_enbi) {}
  virtual vector operator()(const vector& arg);
};
vector EnbiResidualFunc::operator()(const vector& arg) {
  return enbi->getResidual(arg);
}

 
void test_recon() {
  //LinearEnbiproDummy enbi_dummy;
  QuadraticEnbiproDummy enbi_dummy;
  vector x = enbi_dummy.getValues();
  matrix S_x = enbi_dummy.getCovarianceMatrix();
  EnbiResidualFunc f(&enbi_dummy);
  EnbiJacobianFunc J(&enbi_dummy);
  vector v(x.size()), v_exp(5, mischer_teiler::v);;
  v_exp.set(4, -2.5); // note that x(4) is different here, thus v(4) also
  matrix S_v(v.size(), v.size());

  recon(x, S_x, f, J, v, S_v);
  assert(f(x+v).norm2()<1e-6);
  assert_vector_almost_equal_upto(v, v_exp, 1e-3);
}

void recon_test_suite() {
  test_lin_recon_update();
  test_lin_recon();
  test_recon();
}
