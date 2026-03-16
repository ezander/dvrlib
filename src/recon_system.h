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
  void add_covariance_coeff(const char* name1, const char* name2,
                            double cov_coeff);

  int find_var(const std::string& str) const;
  void change_var(const char* name, double val, double confint);
  int get_number_measured() const;
  matrix get_covariance_matrix() const;
  vector get_values() const;
  recon_system updated(const vector& values, const vector& confints) const;
  void print_vars() const;
  void print_constraints(const matrix& F) const;
};

}  // namespace dvrlib

#endif  // DVRLIB_RECON_SYSTEM_H
