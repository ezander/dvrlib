#ifndef __RECON_H__
#define __RECON_H__

namespace dvrlib{

class vector;
class matrix;

template<class argtype, class restype>
class func {
public:
  virtual restype operator()(const argtype& arg) = 0;
};


void lin_cov_update(const matrix& S_x,
		    const matrix& F,
		    matrix& S_v);

void lin_recon(const vector& r,
	       const matrix& S_x,
	       const matrix& F,
	       vector& v);

void lin_recon_update(const vector& r,
		      const matrix& S_x_inv,
		      const matrix& F,
		      const vector& v,
		      vector& dv);

int recon(const vector& x,
	  const matrix& S_x,
	  func<vector, vector>& f,
	  func<vector, matrix>& J,
	  vector& v,
	  matrix& S_v,
	  double eps=1e-6,
	  int maxiter=50);


void extract_confidence(const matrix& S_xnew,
			vector& conf_results);

/**
   Convert a 95% confidence interval into a variance.
*/
double confint2var(double confint);

/**
   Convert a variance into a 95% confidence interval.
*/
double var2confint(double var);

} // namespace dvrlib

#endif // __RECON_H__
