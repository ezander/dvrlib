#include "gsl_wrapper.h"
#include "recon.h"
#include "utils.h"
#include "recon_system.h"

#include <cmath>
#include <string>
#include <vector>
using std::string;

double confint2var(double confint) {
  return pow((confint/1.96),2);
}
double var2confint(double var) {
  return 1.96*sqrt(var);
}

void recon_system::add_var(const char* name, double val, double confint) {
  var v = {name, val, confint};
  vars.push_back(v);
}

void recon_system::add_covariance_coeff(const char* name1, const char* name2, 
					double cov_coeff) {
  extra_cov ec = {name1, name2, cov_coeff};
  extra_covs.push_back(ec);
}

int recon_system::find_var(const string& str) {
  int n = vars.size();
  for( int i=0; i<n; i++){
    if(vars[i].name==str)
      return i;
  }
  return -1;
}
  
int recon_system::get_number_measured() {
  int count = 0, n = vars.size();
  for( int i=0; i<n; i++){
    if(vars[i].confint>=0)
      count++;
  }
  return count;
}

matrix recon_system::get_covariance_matrix() {
  int n = get_number_measured();
  matrix S_x(n, n);
  for(int i=0; i<n; i++){
    S_x.set(i, i, confint2var(vars[i].confint));
  }

  for(unsigned int k=0; k<extra_covs.size(); k++) {
    int i = find_var(extra_covs[k].var1);
    int j = find_var(extra_covs[k].var2);
    double rho = extra_covs[k].cov_coeff;
    double cov_ii = S_x.get(i, i);
    double cov_jj = S_x.get(j, j);
    double cov_ij = rho * sqrt(cov_ii * cov_jj);
    S_x.set(i, j, cov_ij);
    S_x.set(j, i, cov_ij);
  }
  return S_x;
}

vector recon_system::get_values() {
  int n = vars.size();
  vector x(n);
  for( int i=0; i<n; i++){
    x.set(i, vars[i].value);
  }
  return x;
}
  

void recon_system::print_constraints(const matrix& F) {
  for(int i=0; i<F.size1(); i++){
    bool first=true;
    for(int j=0; j<F.size2(); j++) {
      double val = F.get(i,j);
      if(val==0) continue;
      if(val>0) {
	if(!first) std::cout << " + ";
      }
      else {
	if(!first) 
	  std::cout << " - ";
	else
	  std::cout << "-";
	val = -val;
      }
      if(val!=1) std::cout << val << "*";
      std::cout << vars[j].name;
      first = false;
    }
    std::cout << " = 0" << std::endl;
  }
}
