#include "gsl_wrapper.h"
#include "utils.h"

#include <cmath>
#include <cassert>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

namespace dvrlib{

void gsl_err_handler(const char * reason,
		      const char * file,
		      int line,
		      int gsl_errno) {
    gsl_exception ex = {reason, file, line, gsl_errno};
  throw ex;
}



void gsl_enable_exceptions(){
  gsl_set_error_handler(gsl_err_handler);
}


////////////////////////////////////////////////////////////////////////////
// vector class
////////////////////////////////////////////////////////////////////////////
vector::vector() {
  v = 0; // private ctor only for class vector_view
}

vector::vector(int n) {
  v = gsl_vector_alloc(n);
  gsl_vector_set_zero(v);
}

vector::vector(int n, double x) {
  v = gsl_vector_alloc(n);
  gsl_vector_set_all(v, x);
}

vector::vector(int n, const double* x) {
  v = gsl_vector_alloc(n);
  gsl_vector_const_view  src = gsl_vector_const_view_array(x, n);
  gsl_vector_memcpy(v, &src.vector);
}


vector::vector(const vector& src) {
  v = gsl_vector_alloc(src.size());
  gsl_vector_memcpy(v, src.v);
}

vector::~vector() {
  if(v->owner)
    gsl_vector_free(v);
}

gsl_vector* vector::gsl_internal() {
  return v;
}

const gsl_vector* vector::gsl_internal() const {
  return v;
}

int vector::size() const {
  return v->size;
}

void vector::set(int i, double val) {
  gsl_vector_set(v, i, val);
}

double vector::get(int i) const {
  return gsl_vector_get(v, i);
}

vector& vector::operator=(const vector& src) {
  if( &src == this)
    return *this;
  gsl_vector_memcpy(v, src.v);
  return *this;
}

vector& vector::operator+=(const vector& src) {
  gsl_vector_add(v, src.v);
  return *this;
}

vector vector::operator+(const vector& src) const {
  vector vec(*this);
  gsl_vector_add(vec.v , src.v);
  return vec;
}

vector& vector::operator-=(const vector& src) {
  gsl_vector_sub(v, src.v);
  return *this;
}

vector vector::operator-(const vector& src) const {
  vector vec(*this);
  gsl_vector_sub(vec.v , src.v);
  return vec;
}

vector vector::operator-() const {
  vector vec(*this); // -z
  vec = vec * (-1);
  return vec;
}


vector& vector::operator*=(double d) {
  gsl_vector_scale(v, d);
  return *this;
}

vector vector::operator*(double d) const {
  vector vec(*this);
  gsl_vector_scale(vec.v , d);
  return vec;
}

vector_view vector::subvector(int k, int n) {
  vector_view vv(gsl_vector_subvector(v, k, n));
  return vv;
}

const vector_view vector::subvector(int k, int n) const {
  vector_view vv(gsl_vector_subvector(v, k, n));
  return vv;
}

vector operator*(double d, const vector& src){
  return src * d;
}


double vector::norm1() const {
  return gsl_blas_dasum(v);
}

double vector::norm2() const {
  return gsl_blas_dnrm2(v);
}



std::ostream& operator<<(std::ostream& out, const vector& vec){
  out << "[";
  for(int i=0; i<vec.size(); i++) {
    if(i) out << ", ";
    out << vec.get(i);
  }
  out << "]";
  return out;
}

////////////////////////////////////////////////////////////////////////////
// vector_view class
////////////////////////////////////////////////////////////////////////////

vector_view::vector_view(gsl_vector_view _vv) {
  vv = _vv;
  v = &vv.vector;
}

vector_view::vector_view(const vector_view& src) {
  vv = src.vv;
  v = &vv.vector;
}

gsl_vector_view* vector_view::gsl_internal() {
  return &vv;
}

vector_view& vector_view::operator=(const vector& src) {
  if( &src == this)
    return *this;
  gsl_vector_memcpy(v, src.v);
  return *this;
}



////////////////////////////////////////////////////////////////////////////
// matrix class
////////////////////////////////////////////////////////////////////////////


matrix::matrix() {
  m = 0; // private constructor for matrix_view only
}

matrix::matrix(int n1, int n2, bool id, const double* diag) {
  m = gsl_matrix_alloc(n1, n2);
  if(!id)
    gsl_matrix_set_zero(m);
  else
    gsl_matrix_set_identity(m);
  if(diag) {
    for(int i=0; i<n1 && i<n2; i++)
      set(i, i, diag[i]);

  }
}

matrix::matrix(int n1, int n2, const double* x) {
  m = gsl_matrix_alloc(n1, n2);
  gsl_matrix_const_view  src = gsl_matrix_const_view_array(x, n1, n2);
  gsl_matrix_memcpy(m, &src.matrix);
}

matrix::matrix(const matrix& src) {
  m = gsl_matrix_alloc(src.size1(), src.size2());
  gsl_matrix_memcpy(m, src.m);
}

matrix::matrix(gsl_matrix* src) {
    m = gsl_matrix_alloc(src->size1, src->size2);
    gsl_matrix_memcpy(m, src);
}


matrix::~matrix() {
  if(m->owner)
    gsl_matrix_free(m);
}

gsl_matrix* matrix::gsl_internal() {
  return m;
}

const gsl_matrix* matrix::gsl_internal() const {
  return m;
}

int matrix::size1() const {
  return m->size1;
}

int matrix::size2() const {
  return m->size2;
}

void matrix::set(int i, int j, double val) {
  gsl_matrix_set(m, i, j, val);
}

double matrix::get(int i, int j) const {
    return gsl_matrix_get(m, i, j);
}

matrix& matrix::operator=(const matrix& src) {
  if( &src == this)
    return *this;
  gsl_matrix_memcpy(m, src.m);
  return *this;
}

matrix matrix::operator+(const matrix& src) const {
  matrix c(*this);
  assert(c.size1()==src.size1());
  assert(c.size2()==src.size2());

  gsl_matrix_add (c.m , src.m);
  return c;
}

matrix matrix::operator+=(const matrix& src) const {
  gsl_matrix_add (m , src.m);
  return *this;
}

matrix matrix::operator-(const matrix& src) const {
  matrix c(*this);
  assert(c.size1()==src.size1());
  assert(c.size2()==src.size2());

  gsl_matrix_sub (c.m , src.m);
  return c;
}

matrix matrix::operator-=(const matrix& src) const {
  gsl_matrix_sub (m , src.m);
  return *this;
}

matrix matrix::operator-() const {
  matrix c(*this);
  gsl_matrix_scale(c.m, -1.0);
  return c.m;
}


vector matrix::operator*(const vector& src) const {
  vector y(size1());
  assert(size2()==src.size());
  gsl_blas_dgemv(CblasNoTrans, 1.0, m, src.v, 0.0, y.v);
  return y;
}

matrix matrix::operator*(const matrix& src) const {
  matrix c(size1(), src.size2());
  assert(size2()==src.size1());
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, m, src.m, 0.0, c.m);
  return c;
}

matrix matrix::operator*(double d) const {
  matrix c(*this);
  gsl_matrix_scale(c.m, d);
  return c.m;
}


matrix matrix::operator*=(double d) const {
  gsl_matrix_scale(m, d);
  return m;
}




matrix matrix::transpose() const {
  matrix dest(size2(), size1());
  gsl_matrix_transpose_memcpy(dest.m, m);
  return dest;
}

struct permutation {
  gsl_permutation* p;
  permutation(int n) {
    p = gsl_permutation_alloc(n);
  }
  ~permutation() {
    gsl_permutation_free(p);
  }
};

matrix matrix::inverse() const {
  matrix LU(*this), B(*this);
  permutation p(size1());
  int sign;
  gsl_linalg_LU_decomp(LU.m, p.p, &sign);
  gsl_linalg_LU_invert(LU.m, p.p, B.m);
  return B;
}


vector matrix::linsolve(const vector& b) const {
  matrix LU(*this);
  vector x(b.size());
  permutation p(size1());
  int sign;
  gsl_linalg_LU_decomp(LU.m, p.p, &sign);
  gsl_linalg_LU_solve(LU.m, p.p, b.gsl_internal(), x.v);
  return x;
}

void matrix::svd(matrix& U, matrix& V, vector& s) const {
  int M = size1();
  int N = size2();
  assert(M>=N);
  U = *this;
  vector work(N);
  gsl_linalg_SV_decomp(U.m, V.m, s.v, work.v);
}

vector matrix::svd() const {
  int M = size1();
  int N = size2();
  if(M<N) 
    return transpose().svd();
  vector s(N);
  matrix U(M, N);
  matrix V(N, N);
  svd(U, V, s);
  return s;
}


matrix_view matrix::submatrix(int k1, int k2, int n1, int n2) {
  assert(k1+n1<=size1());
  assert(k2+n2<=size2());
  matrix_view mv(gsl_matrix_submatrix(m, k1, k2, n1, n2));
  return mv;
}

const matrix_view matrix::submatrix(int k1, int k2, int n1, int n2) const {
  assert(k1+n1<=size1());
  assert(k2+n2<=size2());
  matrix_view mv(gsl_matrix_submatrix(m, k1, k2, n1, n2));
  return mv;
}

std::ostream& operator<<(std::ostream& out, const matrix& mat){
  for(int i=0; i<mat.size1(); i++) {
    if(i) 
      out << "]" << std::endl << " [";
    else
      out << "[[";
    for(int j=0; j<mat.size2(); j++) {
      if(j) out << ", ";
      out << mat.get(i,j);
    }
  }
  out << "]]";
  return out;
}



////////////////////////////////////////////////////////////////////////////
// matrix_view class
////////////////////////////////////////////////////////////////////////////

matrix_view::matrix_view(gsl_matrix_view _mv) {
  mv = _mv;
  m = &mv.matrix;
}

matrix_view::matrix_view(const matrix_view& src) {
  mv = src.mv;
  m = &mv.matrix;
}

gsl_matrix_view* matrix_view::gsl_internal() {
  return &mv;
}


matrix_view& matrix_view::operator=(const matrix& src) {
  if( &src == this)
    return *this;
  gsl_matrix_memcpy(m, src.m);
  return *this;
}


matrix operator*(double d, const matrix& src){
    gsl_matrix* dest = gsl_matrix_alloc(src.size1(), src.size2());
    gsl_matrix_memcpy(dest, src.gsl_internal());
    gsl_matrix_scale(dest, d);
    return dest;
}


////////////////////////////////////////////////////////////////////////////
// gsl helpers
////////////////////////////////////////////////////////////////////////////

#define PRINT_START(cls)\
  out << cls << " {" << std::endl
#define PRINT_FIELD(s,field)\
  out << "  " << #field << ": " << s.field << std::endl
#define PRINT_END()\
  out << "}" << std::endl

std::ostream& operator<<(std::ostream& out, const gsl_vector& v) {
  PRINT_START("gls_vector");
  PRINT_FIELD(v, size);
  PRINT_FIELD(v, stride);
  PRINT_FIELD(v, data);
  PRINT_FIELD(v, block);
  PRINT_FIELD(v, owner);
  PRINT_END();
  return out;
}

std::ostream& operator<<(std::ostream& out, const gsl_matrix& m) {
  PRINT_START("gls_matrix");
  PRINT_FIELD(m, size1);
  PRINT_FIELD(m, size2);
  PRINT_FIELD(m, tda);
  PRINT_FIELD(m, data);
  PRINT_FIELD(m, block);
  PRINT_FIELD(m, owner);
  PRINT_END();
  return out;
}
std::ostream& operator<<(std::ostream& out, const gsl_matrix_view& mv) {
  PRINT_START("gls_matrix_view");
  PRINT_FIELD(mv.matrix, size1);
  PRINT_FIELD(mv.matrix, size2);
  PRINT_FIELD(mv.matrix, tda);
  PRINT_FIELD(mv.matrix, data);
  PRINT_FIELD(mv.matrix, block);
  PRINT_FIELD(mv.matrix, owner);
  PRINT_END();
  return out;
}

}
