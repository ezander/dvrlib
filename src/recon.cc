#include <cassert>
#include "recon.h"
#include "gsl_wrapper.h"
#include "utils.h"

void lin_cov_update_Streit(const matrix& S_x, const matrix& F, matrix& S_v) {
  int M = S_x.size1();
  int N = F.size2() - M;
  int K = F.size1();

  assert(S_x.size2() == M);
  assert(N>=0);

  matrix Fs = F.submatrix(0, 0, K, M);

  matrix S_x_F_T = S_x * Fs.transpose();
  S_v = S_x_F_T * (Fs * S_x_F_T).inverse() * S_x_F_T.transpose();
}

/** 
 * Compute the update of the covariance matrix. The update is computed by
 * the following formula
 * \f[ S_v = G S_F G^T \f]
 * where 
 * \f[ G = P_m * Z^{-1} * P_l \f]
 * and
 * \f[ S_F = F_m * S_x * F_m. \f]
 * \f$P_m\f$ and \f$P_l\f$ indicate projection matrices on the space of 
 * measured variables and on the space of Lagrange multipliers, respectively.
 * 
 * @param[in] S_x The input covariance matrix
 * @param[in] F The matrix of constraints
 * @param[out] S_v The covariance of the updates
 */
void lin_cov_update_Zander(const matrix& S_x, const matrix& F, matrix& S_v) {
  int M = S_x.size1();
  int N = F.size2() - M;
  int K = F.size1();
  int L = M + N + K;
  
  assert(S_x.size2() == M);
  assert(N>=0);

  matrix P_m(M, L);
  matrix I_m(M, M, true);
  P_m.submatrix(0, 0, M, M) = I_m;

  matrix P_l(K, L);
  matrix I_l(K, K, true);
  P_l.submatrix(0, M+N, K, K) = I_l;

  matrix Z(L, L);
  Z.submatrix(0, 0, M, M) = S_x.inverse();
  Z.submatrix(0, M+N, M+N, K) = F.transpose();
  Z.submatrix(M+N, 0, K, M+N) = F;
  matrix Z_inv = Z.inverse();

  matrix F_m = F.submatrix(0, 0, K, M);

  matrix G = P_m * Z_inv * P_l.transpose();
  matrix S_F = F_m * S_x * F_m.transpose();

  S_v = G * S_F * G.transpose();
}



void lin_cov_update(const matrix& S_x, const matrix& F, matrix& S_v) {
  //lin_cov_update_Streit(S_x, F, S_v);
  lin_cov_update_Zander(S_x, F, S_v);
}


void lin_recon(const vector& r,
	       const matrix& S_x,
	       const matrix& F,
	       vector& v) {

  int M = S_x.size1();
  int N = F.size2() - M;
  int K = F.size1();

  assert(S_x.size2() == M);
  assert(r.size() == K);
  assert(v.size() == M+N);

  matrix Z(M+N+K, M+N+K);
  Z.submatrix(0, 0, M, M) = S_x.inverse();
  Z.submatrix(0, M+N, M+N, K) = F.transpose();
  Z.submatrix(M+N, 0, K, M+N) = F;

  vector g(M+N+K);
  g.subvector(M+N, K) = -1.0 * r;

  vector z = Z.linsolve(g);
  v = z.subvector(0, M+N);

  matrix S_v(S_x);
  lin_cov_update(S_x, F, S_v);
}


void lin_recon_update(const vector& r,
		      const matrix& S_x_inv,
		      const matrix& F,
		      const vector& v,
		      vector& dv) {

  int M = S_x_inv.size1();
  int N = v.size() - M;
  int K = r.size();

  assert(S_x_inv.size2() == M);
  assert(F.size1() == K);
  assert(F.size2() == M+N);
  assert(dv.size() == M+N);

  matrix Z(M+N+K, M+N+K);
  Z.submatrix(0, 0, M, M) = S_x_inv;
  Z.submatrix(0, M+N, M+N, K) = F.transpose();
  Z.submatrix(M+N, 0, K, M+N) = F;

  vector g(M+N+K);
  g.subvector(0, M) = -1.0 * (S_x_inv * v.subvector(0, M));
  g.subvector(M+N, K) = -1.0 * r;
  
  vector z = Z.linsolve(g);
  dv = z.subvector(0, M+N);
}


int recon(const vector& x,
	   const matrix& S_x, 
	   func<vector, vector>& f,
	   func<vector, matrix>& J,
	   vector& v, 
	   matrix& S_v,
	   double eps,
	   int maxiter) {

  matrix S_x_inv = S_x.inverse();

  for(int iter=0; iter<maxiter+1; iter++) {
    if( iter==maxiter) {
      std::cout << "recon did not converge" << std::endl;
      return 1;
    }

    vector r=f(x + v);
    // exit is residual is small enough
    if(r.norm2()<eps) {
      break;
    }

    matrix F = J(x + v);
    vector dv(0*v);
    
    lin_recon_update(r, S_x_inv, F, v, dv);

    v = v + dv;
    // update to v is smaller than eps (should be a different eps 
    // than that above...)
    if(dv.norm2()<eps) {
      break;
    }
  }
  
  return 0;
}
