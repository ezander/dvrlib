#ifndef DVRLIB_GSL_WRAPPER_H
#define DVRLIB_GSL_WRAPPER_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <iostream>
#include "dvr_assert.h"

namespace dvrlib {

/** Exception type thrown by GSL when \c gsl_enable_exceptions() is active. */
struct gsl_exception {
  const char* reason;  ///< Human-readable error description
  const char* file;    ///< Source file where the error occurred
  int line;            ///< Line number where the error occurred
  int gsl_errno;       ///< GSL error code (see \c gsl_errno.h)
};

/** Install a GSL error handler that throws \c gsl_exception instead of aborting. */
void gsl_enable_exceptions();

/**
 * Switch between standard and VDI 2048 pretty-print format for vectors and
 * matrices written to \c std::ostream.
 *
 * @param[in] enable \c true selects VDI format, \c false selects standard format.
 */
void set_vdi_print_format(bool enable);

class vector_view;

/**
 * RAII wrapper around a heap-allocated \c gsl_vector.
 *
 * Owns the underlying GSL allocation and frees it on destruction.
 * Arithmetic operators follow GSL semantics (element-wise).
 */
class vector {
  gsl_vector* v;

public:
  /** Allocate an uninitialised vector of length \c n. */
  vector(int n);
  /** Allocate a vector of length \c n with all elements set to \c x. */
  vector(int n, double x);
  /** Allocate a vector of length \c n and copy values from the array \c x. */
  vector(int n, const double* x);
  /** Copy constructor. */
  vector(const vector& src);
  /** Destructor — frees the underlying GSL allocation. */
  ~vector();

  /** Return a pointer to the underlying \c gsl_vector (mutable). */
  gsl_vector* gsl_internal();
  /** Return a pointer to the underlying \c gsl_vector (read-only). */
  const gsl_vector* gsl_internal() const;

  /** Return the number of elements. */
  int size() const;
  /** Set element \c i to \c val. */
  void set(int i, double val);
  /** Return element \c i. */
  double get(int i) const;
  /** Return element \c i (read-only subscript). */
  double operator[](int i) const;

  /** Copy-assign from another vector. */
  vector& operator=(const vector& src);
  /** Copy-assign from a view. */
  vector& operator=(const vector_view& src);
  vector& operator+=(const vector& src);
  vector operator+(const vector& src) const;
  vector& operator-=(const vector& src);
  vector operator-(const vector& src) const;
  /** Unary negation — returns \c -(*this). */
  vector operator-() const;

  /** Scale all elements by \c d in place. */
  vector& operator*=(double d);
  /** Return a copy scaled by \c d. */
  vector operator*(double d) const;
  /** Return the L1 norm (sum of absolute values). */
  double norm1() const;
  /** Return the L2 (Euclidean) norm. */
  double norm2() const;

  /**
   * Return a view of \c n elements starting at index \c k.
   *
   * @param[in] k Start index (0-based).
   * @param[in] n Number of elements in the sub-vector.
   */
  vector_view subvector(int k, int n);
  /** @copydoc vector::subvector(int,int) */
  const vector_view subvector(int k, int n) const;
  friend class matrix;
};

/** Scalar multiplication \c d*src (commutative with \c vector::operator*). */
vector operator*(double d, const vector& src);
/** Write a vector to an output stream. Format depends on \c set_vdi_print_format(). */
std::ostream& operator<<(std::ostream& out, const vector& vec);

/**
 * Non-owning view into a contiguous slice of a \c vector or \c matrix row/column.
 *
 * Wraps a \c gsl_vector_view. The referenced data is owned by the original
 * \c vector or \c matrix — the view must not outlive it.
 * Assigning to a view modifies the underlying data in place.
 * Use \c operator vector() to obtain an independent copy.
 */
class vector_view {
  gsl_vector_view vv;
  vector_view(const gsl_vector_view& vv);

public:
  /** Copy constructor (shallow — both views reference the same data). */
  vector_view(const vector_view& src);

  /** Return the number of elements. */
  int size() const;
  /** Return element \c i. */
  double get(int i) const;
  /** Set element \c i to \c val (modifies the underlying data). */
  void set(int i, double val);
  /** Return a read-only pointer to the underlying \c gsl_vector. */
  const gsl_vector* gsl_internal() const;

  /** Copy values from \c src into the viewed memory. */
  vector_view& operator=(const vector& src);
  /** Copy values from another view into the viewed memory. */
  vector_view& operator=(const vector_view& src);
  /** Create an independent \c vector copy of this view. */
  operator vector() const;

  friend class vector;
  friend class matrix;
};

class matrix_view;

/**
 * RAII wrapper around a heap-allocated \c gsl_matrix.
 *
 * Owns the underlying GSL allocation and frees it on destruction.
 * Row access via \c operator[] returns a \c vector_view into the matrix data.
 */
class matrix {
  gsl_matrix* m;

public:
  /**
   * Allocate an \c n1 × \c n2 matrix.
   *
   * @param[in] n1   Number of rows.
   * @param[in] n2   Number of columns.
   * @param[in] id   If \c true, initialise as identity (requires \c n1 == \c n2
   *                 when \c diag is null).
   * @param[in] diag Optional array of \c n1 diagonal values; if non-null
   *                 overrides \c id.
   */
  matrix(int n1, int n2, bool id = false, const double* diag = nullptr);
  /** Allocate an \c n1 × \c n2 matrix and copy values from the row-major array \c x. */
  matrix(int n1, int n2, const double* x);
  /**
   * Allocate an \c n1 × \c n2 matrix from a 2-D C array.
   *
   * The template parameter \c n must match \c n2.
   */
  template <int n>
  matrix(int n1, int n2, const double (*x)[n]);

  /** Copy constructor. */
  matrix(const matrix& src);
  /** Construct from an existing (non-owned) \c gsl_matrix pointer. */
  matrix(const gsl_matrix* src);
  /** Destructor — frees the underlying GSL allocation. */
  ~matrix();

  /** Return a pointer to the underlying \c gsl_matrix (mutable). */
  gsl_matrix* gsl_internal();
  /** Return a pointer to the underlying \c gsl_matrix (read-only). */
  const gsl_matrix* gsl_internal() const;

  /** Return the number of rows. */
  int size1() const;
  /** Return the number of columns. */
  int size2() const;
  /** Set element (\c i, \c j) to \c val. */
  void set(int i, int j, double val);
  /** Return element (\c i, \c j). */
  double get(int i, int j) const;
  /** Return a mutable view of row \c i. */
  vector_view operator[](int i);
  /** Return a read-only view of row \c i. */
  const vector_view operator[](int i) const;

  /** Copy-assign from another matrix. */
  matrix& operator=(const matrix& src);
  /** Copy-assign from a view. */
  matrix& operator=(const matrix_view& src);

  matrix operator+(const matrix& src) const;
  matrix& operator+=(const matrix& src);
  matrix operator-(const matrix& src) const;
  matrix& operator-=(const matrix& src);
  /** Unary negation — returns \c -(*this). */
  matrix operator-() const;

  /** Matrix–vector product. */
  vector operator*(const vector& src) const;
  /** Matrix–matrix product. */
  matrix operator*(const matrix& src) const;
  /** Return a copy scaled by \c d. */
  matrix operator*(double d) const;
  /** Scale all elements by \c d in place. */
  matrix& operator*=(double d);

  /** Return the transpose. */
  matrix transpose() const;
  /** Return the inverse (via LU decomposition). */
  matrix inverse() const;
  /**
   * Solve the linear system \c A*x = b and return \c x.
   *
   * @param[in] b Right-hand side vector (length \c size1()).
   */
  vector linsolve(const vector& b) const;
  /**
   * Compute the full singular value decomposition \c A = U * diag(s) * V^T.
   *
   * @param[out] U Left singular vectors (same size as \c *this).
   * @param[out] V Right singular vectors (\c size2() × \c size2()).
   * @param[out] s Singular values (length \c size2()).
   */
  void svd(matrix& U, matrix& V, vector& s) const;
  /** Compute and return only the singular values. */
  vector svd() const;

  /**
   * Return a view of the \c n1 × \c n2 sub-matrix starting at (\c k1, \c k2).
   *
   * @param[in] k1 Starting row index (0-based).
   * @param[in] k2 Starting column index (0-based).
   * @param[in] n1 Number of rows in the sub-matrix.
   * @param[in] n2 Number of columns in the sub-matrix.
   */
  matrix_view submatrix(int k1, int k2, int n1, int n2);
  /** @copydoc matrix::submatrix(int,int,int,int) */
  const matrix_view submatrix(int k1, int k2, int n1, int n2) const;
  friend class matrix_view;
};

template <int n>
matrix::matrix(int n1, int n2, const double (*x)[n]) {
  dvr_assert(n == n2);
  m = gsl_matrix_alloc(n1, n2);
  gsl_matrix_const_view src = gsl_matrix_const_view_array(x[0], n1, n2);
  gsl_matrix_memcpy(m, &src.matrix);
}

/** Scalar multiplication \c d*src (commutative with \c matrix::operator*). */
matrix operator*(double d, const matrix& src);

/** Write a matrix to an output stream. Format depends on \c set_vdi_print_format(). */
std::ostream& operator<<(std::ostream& out, const matrix& mat);

/**
 * Non-owning view into a rectangular sub-region of a \c matrix.
 *
 * Wraps a \c gsl_matrix_view. The referenced data is owned by the original
 * \c matrix — the view must not outlive it.
 * Assigning to a view modifies the underlying data in place.
 * Use \c operator matrix() to obtain an independent copy.
 */
class matrix_view {
  gsl_matrix_view mv;
  matrix_view(const gsl_matrix_view& mv);

public:
  /** Copy constructor (shallow — both views reference the same data). */
  matrix_view(const matrix_view& src);

  /** Return the number of rows. */
  int size1() const;
  /** Return the number of columns. */
  int size2() const;
  /** Return element (\c i, \c j). */
  double get(int i, int j) const;
  /** Set element (\c i, \c j) to \c val (modifies the underlying data). */
  void set(int i, int j, double val);
  /** Return a read-only pointer to the underlying \c gsl_matrix. */
  const gsl_matrix* gsl_internal() const;

  /** Copy values from \c src into the viewed memory. */
  matrix_view& operator=(const matrix& src);
  /** Copy values from another view into the viewed memory. */
  matrix_view& operator=(const matrix_view& src);
  /** Create an independent \c matrix copy of this view. */
  operator matrix() const;

  friend class matrix;
};

std::ostream& operator<<(std::ostream& out, const gsl_vector& v);
std::ostream& operator<<(std::ostream& out, const gsl_matrix& m);
std::ostream& operator<<(std::ostream& out, const gsl_matrix_view& mv);

}  // namespace dvrlib

#endif  // DVRLIB_GSL_WRAPPER_H
