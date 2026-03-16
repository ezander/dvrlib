#ifndef DVRLIB_RECON_H
#define DVRLIB_RECON_H

namespace dvrlib {

class vector;
class matrix;

/**
 * Abstract function object mapping \c argtype to \c restype.
 *
 * Derive from this to supply residual functions and Jacobians to \c recon().
 */
template <class argtype, class restype>
class func {
public:
  virtual ~func() = default;
  /** Evaluate the function at \c arg. */
  virtual restype operator()(const argtype& arg) const = 0;
};

/**
 * Propagate measurement covariance through a linear constraint matrix.
 *
 * Computes \c S_v = F * S_x * F^T, the covariance of the reconciled
 * values given the measurement covariance \c S_x and the constraint
 * matrix \c F.
 *
 * @param[in]  S_x Measurement covariance matrix.
 * @param[in]  F   Linear constraint (Jacobian) matrix.
 * @param[out] S_v Resulting covariance of the reconciled values.
 */
void lin_cov_update(const matrix& S_x,
                    const matrix& F,
                    matrix& S_v);

/**
 * Solve the linear reconciliation problem.
 *
 * Computes the reconciled values \c v that satisfy the linear constraints
 * encoded in \c F, given raw measurements \c r and measurement covariance
 * \c S_x.
 *
 * @param[in]  r   Raw measurement vector.
 * @param[in]  S_x Measurement covariance matrix.
 * @param[in]  F   Linear constraint matrix.
 * @param[out] v   Reconciled values.
 */
void lin_recon(const vector& r,
               const matrix& S_x,
               const matrix& F,
               vector& v);

/**
 * Compute the reconciliation correction step for non-linear iteration.
 *
 * Given the current iterate \c v, computes the Newton-like correction
 * \c dv to reduce the residual of the constraint equations.
 *
 * @param[in]  r       Raw measurement vector.
 * @param[in]  S_x_inv Inverse of the measurement covariance matrix.
 * @param[in]  F       Constraint Jacobian evaluated at the current iterate.
 * @param[in]  v       Current reconciled values.
 * @param[out] dv      Correction to apply: next iterate is \c v + dv.
 */
void lin_recon_update(const vector& r,
                      const matrix& S_x_inv,
                      const matrix& F,
                      const vector& v,
                      vector& dv);

/**
 * Solve the non-linear reconciliation problem by iterative linearisation.
 *
 * Iterates a Newton-type update until the correction \c dv is smaller than
 * \c eps (in the L2 norm) or \c maxiter iterations are reached.
 *
 * @param[in]  x       Raw measurement vector.
 * @param[in]  S_x     Measurement covariance matrix.
 * @param[in]  f       Constraint residual function \c f(v) = 0.
 * @param[in]  J       Jacobian of \c f with respect to \c v.
 * @param[out] v       Reconciled values.
 * @param[out] S_v     Covariance of the reconciled values.
 * @param[in]  eps     Convergence tolerance on the correction norm (default 1e-6).
 * @param[in]  maxiter Maximum number of iterations (default 50).
 * @return Number of iterations performed, or -1 if the method did not converge.
 */
int recon(const vector& x,
          const matrix& S_x,
          const func<vector, vector>& f,
          const func<vector, matrix>& J,
          vector& v,
          matrix& S_v,
          double eps = 1e-6,
          int maxiter = 50);

/**
 * Extract 95% confidence intervals from a covariance matrix diagonal.
 *
 * @param[in]  S_xnew       Covariance matrix of the reconciled values.
 * @param[out] conf_results Confidence intervals (one per variable).
 */
void extract_confidence(const matrix& S_xnew,
                        vector& conf_results);

/**
 * Convert a 95% confidence interval into a variance.
 *
 * Uses the relation \c confint = z_{0.975} * sqrt(var) where
 * \c z_{0.975} = Φ⁻¹(0.975) ≈ 1.96.
 */
double confint2var(double confint);

/**
 * Convert a variance into a 95% confidence interval.
 *
 * Uses the relation \c confint = z_{0.975} * sqrt(var) where
 * \c z_{0.975} = Φ⁻¹(0.975) ≈ 1.96.
 */
double var2confint(double var);

}  // namespace dvrlib

#endif  // DVRLIB_RECON_H
