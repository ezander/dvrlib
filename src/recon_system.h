#ifndef DVRLIB_RECON_SYSTEM_H
#define DVRLIB_RECON_SYSTEM_H

#include <vector>
#include <string>
#include "gsl_wrapper.h"

namespace dvrlib {

/**
   Encapsulates a system to be reconciled, including variables,
   confidence intervals, covariance coefficient, etc. Constraints still
   need to be implemented. Also, fixed variables should be treated
   here.
 */
class recon_system {
  struct var {
    std::string name;
    double value;
    double confint;
  };
  struct extra_cov {
    std::string var1;
    std::string var2;
    double cov_coeff;
  };

  std::vector<var> vars;
  std::vector<extra_cov> extra_covs;

public:
  /**
   * Add a variable to be reconciled to the system. If \c confint is
   * negative, it is assumed that this variable is free (i.e., no
   * measurements). \todo If \c confint is zero it is assumed that this
   * variable is fixed.
   *
   * @param[in] name The name of the variable
   * @param[in] val The measured, estimated of fixed value of the variable
   * @param[in] confint The confidence interval
   */
  void add_var(const char* name, double val, double confint);

  /**
   * Register a covariance coefficient between two variables.
   *
   * The coefficient is used when building the covariance matrix via
   * \c get_covariance_matrix().
   *
   * @param[in] name1     Name of the first variable.
   * @param[in] name2     Name of the second variable.
   * @param[in] cov_coeff Covariance coefficient (fraction of the geometric
   *                      mean of the two variances).
   */
  void add_covariance_coeff(const char* name1, const char* name2,
                            double cov_coeff);

  /**
   * Return the index of a variable by name, or -1 if not found.
   *
   * @param[in] str Variable name to search for.
   */
  int find_var(const std::string& str) const;

  /**
   * Update the measured value and confidence interval of an existing variable.
   *
   * @param[in] name    Name of the variable to update.
   * @param[in] val     New measured value.
   * @param[in] confint New confidence interval.
   */
  void change_var(const char* name, double val, double confint);

  /** Return the number of variables that have measurements (non-negative \c confint). */
  int get_number_measured() const;

  /**
   * Build and return the full covariance matrix from the registered variables
   * and covariance coefficients.
   */
  matrix get_covariance_matrix() const;

  /** Return a vector of all variable values in registration order. */
  vector get_values() const;

  /**
   * Return a new \c recon_system with updated values and confidence intervals.
   *
   * Keeps all variable names and covariance coefficients; replaces values and
   * confidence intervals from the supplied vectors.
   *
   * @param[in] values   New variable values (one per variable, in order).
   * @param[in] confints New confidence intervals (one per variable, in order).
   */
  recon_system updated(const vector& values, const vector& confints) const;

  /** Print all variables with their values and confidence intervals to stdout. */
  void print_vars() const;

  /**
   * Print the constraint matrix \c F together with variable names to stdout.
   *
   * @param[in] F Constraint (Jacobian) matrix whose columns correspond to
   *              the variables in this system.
   */
  void print_constraints(const matrix& F) const;
};

}  // namespace dvrlib

#endif  // DVRLIB_RECON_SYSTEM_H
