#include "gsl_wrapper.h"
#include "recon.h"
#include "recon_system.h"
#include "dvr_assert.h"

#include <cstdio>
#include <cmath>
#include <string>

using std::string;

namespace dvrlib {

void recon_system::add_var(const char* name, double val, double confint) {
  var v = {name, val, confint};
  vars.push_back(v);
}

void recon_system::add_covariance_coeff(const char* name1, const char* name2,
                                        double cov_coeff) {
  extra_cov ec = {name1, name2, cov_coeff};
  extra_covs.push_back(ec);
}

int recon_system::find_var(const string& str) const {
  int n = vars.size();
  for(int i = 0; i < n; i++) {
    if(vars[i].name == str)
      return i;
  }
  return -1;
}

void recon_system::change_var(const char* name, double val, double confint) {
  int i = find_var(name);
  dvr_assert(i >= 0);
  vars[i].value = val;
  vars[i].confint = confint;
}

int recon_system::get_number_measured() const {
  int count = 0;
  int n = vars.size();
  for(int i = 0; i < n; i++) {
    if(vars[i].confint >= 0)
      count++;
  }
  return count;
}

matrix recon_system::get_covariance_matrix() const {
  int n = get_number_measured();
  // measured variables must precede free variables
  for(int i = 0; i < n; i++) {
    dvr_assert(vars[i].confint >= 0);
  }
  matrix S_x(n, n);
  for(int i = 0; i < n; i++) {
    S_x.set(i, i, confint2var(vars[i].confint));
  }

  for(unsigned int k = 0; k < extra_covs.size(); k++) {
    int i = find_var(extra_covs[k].var1);
    int j = find_var(extra_covs[k].var2);
    dvr_assert(i >= 0);
    dvr_assert(j >= 0);
    dvr_assert(i < n);
    dvr_assert(j < n);
    double rho = extra_covs[k].cov_coeff;
    double cov_ii = S_x.get(i, i);
    double cov_jj = S_x.get(j, j);
    double cov_ij = rho * sqrt(cov_ii * cov_jj);
    S_x.set(i, j, cov_ij);
    S_x.set(j, i, cov_ij);
  }
  return S_x;
}

vector recon_system::get_values() const {
  int n = vars.size();
  vector x(n);
  for(int i = 0; i < n; i++) {
    x.set(i, vars[i].value);
  }
  return x;
}

recon_system recon_system::updated(const vector& values,
                                   const vector& confints) const {
  dvr_assert(values.size() == (int)vars.size());
  dvr_assert(confints.size() == get_number_measured());
  recon_system result(*this);
  for(int i = 0; i < (int)vars.size(); i++) {
    result.vars[i].value = values.get(i);
  }
  for(int i = 0; i < (int)confints.size(); i++) {
    result.vars[i].confint = confints.get(i);
  }
  return result;
}

void recon_system::print_vars() const {
  printf("%-12s%10s%10s\n", "name", "value", "confint");
  for(unsigned int i = 0; i < vars.size(); i++) {
    printf("%-12s%10.3f", vars[i].name.c_str(), vars[i].value);
    if(vars[i].confint < 0)
      printf("    (free)");
    else
      printf("%10.3f", vars[i].confint);
    printf("\n");
  }
}

void recon_system::print_constraints(const matrix& F) const {
  for(int i = 0; i < F.size1(); i++) {
    bool first = true;
    for(int j = 0; j < F.size2(); j++) {
      double val = F.get(i, j);
      if(val == 0)
        continue;
      if(val > 0) {
        if(!first)
          printf(" + ");
      } else {
        if(!first)
          printf(" - ");
        else
          printf("-");
        val = -val;
      }
      if(val != 1)
        printf("%g*", val);
      printf("%s", vars[j].name.c_str());
      first = false;
    }
    printf(" = 0\n");
  }
}

}  // namespace dvrlib
