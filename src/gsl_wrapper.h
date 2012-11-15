#ifndef __GSL_WRAPPER_H__
#define __GSL_WRAPPER_H__

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <iostream>
#include <cassert>

namespace dvrlib{

struct gsl_exception {
  const char* reason;
  const char* file;
  int line;
  int gsl_errno;
};


void gsl_enable_exceptions();


class vector_view;

class vector {
  gsl_vector* v;
  vector();
public:
  vector(int n);
  template<int n> vector();
  vector(int n, double x);
  vector(int n, const double* x);
  vector(const vector& src);
  ~vector();

  gsl_vector* gsl_internal();
  const gsl_vector* gsl_internal() const;

  int size() const;
  void set(int i, double val);
  double get(int i) const;
  vector& operator=(const vector& src);
  vector& operator+=(const vector& src);
  vector operator+(const vector& src) const;
  vector& operator-=(const vector& src);
  vector operator-(const vector& src) const;
  vector& operator*=(double d);
  vector operator*(double d) const;
  double norm1() const;
  double norm2() const;

  vector_view subvector(int k, int n);
  const vector_view subvector(int k, int n) const;
  friend class vector_view;
  friend class matrix;
};

vector operator*(double d, const vector& src);
std::ostream& operator<<(std::ostream& out, const vector& vec);

class vector_view: public vector {
  gsl_vector_view vv;
  vector_view(gsl_vector_view vv);

public:
  vector_view(const vector_view& src);
  gsl_vector_view* gsl_internal();
  vector_view& operator=(const vector& src);

  friend class vector;
};



class matrix_view;

class matrix {
  gsl_matrix* m;
  matrix();
public:
  matrix(int n1, int n2, bool id=false, const double* diag=0);
  matrix(int n1, int n2, const double* x);
  template<int n>
    matrix(int n1, int n2, const double (*x)[n]);
  
  matrix(const matrix& src);
  matrix(gsl_matrix* src);
  ~matrix();

  gsl_matrix* gsl_internal();
  const gsl_matrix* gsl_internal() const;

  int size1() const;
  int size2() const;
  void set(int i, int j, double val);
  double get(int i, int j) const;
  matrix& operator=(const matrix& src);

  matrix operator+(const matrix& src) const;
  matrix operator+=(const matrix& src) const;
  matrix operator-(const matrix& src) const;
  matrix operator-=(const matrix& src) const;

  vector operator*(const vector& src) const;
  matrix operator*(const matrix& src) const;
 
  matrix transpose() const;
  matrix inverse() const;
  vector linsolve(const vector& b) const;
  void svd(matrix& U, matrix& V, vector& s) const;
  vector svd() const;

  matrix_view submatrix(int k1, int k2, int n1, int n2);
  const matrix_view submatrix(int k1, int k2, int n1, int n2) const;
  friend class matrix_view;
};

template<int n>
matrix::matrix(int n1, int n2, const double (*x)[n]) {
  assert(n==n2);
  m = gsl_matrix_alloc(n1, n2);
  gsl_matrix_const_view  src = gsl_matrix_const_view_array(x[0], n1, n2);
  gsl_matrix_memcpy(m, &src.matrix);
}

std::ostream& operator<<(std::ostream& out, const matrix& vec);

class matrix_view: public matrix {
  gsl_matrix_view mv;
  matrix_view(gsl_matrix_view mv);

public:
  matrix_view(const matrix_view& src);
  gsl_matrix_view* gsl_internal();
  matrix_view& operator=(const matrix& src);

  friend class matrix;
};


std::ostream& operator<<(std::ostream& out, const gsl_vector& v);
std::ostream& operator<<(std::ostream& out, const gsl_matrix& m);
std::ostream& operator<<(std::ostream& out, const gsl_matrix_view& mv);

}


#endif // __GSL_WRAPPER_H__
