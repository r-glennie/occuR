// Author: Richard Glennie
// Date: Feb 2020 
// Occupancy model 

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  // FLAG
  DATA_INTEGER(flag); // flag for computing normalising constant
  
  // DATA 
  DATA_INTEGER(nsites); // number of sites 
  DATA_VECTOR(nocc); // number of occasions 
  DATA_VECTOR(y); // record in order of site then occasion then visit 
  DATA_VECTOR(totsite); // total number of detections per site x occasion 
  DATA_VECTOR(nvisit); // number of visits per site x occasion
  DATA_MATRIX(X_psi); // occupancy fixed effects design matrix
  DATA_MATRIX(Z_psi); // occupancy random effects design matrix 
  DATA_SPARSE_MATRIX(S_psi); // occupancy smoothing matrix 
  DATA_IVECTOR(S_psi_n); // detection dimension of smoothing matrix
  DATA_MATRIX(X_p); // detection fixed effects design matrix
  DATA_MATRIX(Z_p); // detection random effects design matrix 
  DATA_SPARSE_MATRIX(S_p); // detection smoothing matrix 
  DATA_IVECTOR(S_p_n); // detection dimension of smoothing matrix
  
  // PARAMETERS 
  PARAMETER_VECTOR(beta_psi); // fixed effects for occupancy 
  PARAMETER_VECTOR(beta_p); // fixed effects for detection 
  PARAMETER_VECTOR(z_psi); // random effects for occupancy 
  PARAMETER_VECTOR(z_p); // random effects for detection 
  PARAMETER_VECTOR(log_lambda_psi); // log smoothing parameter for occupancy 
  PARAMETER_VECTOR(log_lambda_p); // log smoothing parameter for detection 
  vector<Type> lambda_psi = exp(log_lambda_psi); // smoothing parameter for occupancy 
  vector<Type> lambda_p = exp(log_lambda_p); // smoothing parameter  for detection 
  
  // NEGATIVE-LOG-LIKELIHOOD 
  Type nll = 0; 
  
  // SMOOTHING PENALTIES
  int s0 = 0; 
  // occupancy
  if (S_psi_n(0) > 0) {
    for (int i = 0; i < S_psi_n.size(); ++i) {
      int sn = S_psi_n(i); 
      SparseMatrix<Type> Si = S_psi.block(s0, s0, sn, sn); 
      vector<Type> z_psi_i = z_psi.segment(s0, sn); 
      nll -= Type(0.5) * sn * log_lambda_psi(i) - Type(0.5) * lambda_psi(i) * GMRF(Si, false).Quadform(z_psi_i); 
      s0 += sn; 
    }
  }
  // detection
  if (S_p_n(0) > 0) {
    s0 = 0; 
    for (int i = 0; i < S_p_n.size(); ++i) {
      int sn = S_p_n(i); 
      SparseMatrix<Type> Si = S_p.block(s0, s0, sn, sn); 
      vector<Type> z_p_i = z_p.segment(s0, sn); 
      nll -= Type(0.5) * sn * log_lambda_p(i) - Type(0.5) * lambda_p(i) * GMRF(Si, false).Quadform(z_p_i); 
      s0 += sn; 
    }
  }
  // return un-normalized density for flag
  if (flag == 0) return nll; 
  
  // LINEAR PREDICTORS 
  vector<Type> logit_psi = X_psi * beta_psi; 
  if (S_psi_n(0) > 0) logit_psi += Z_psi * z_psi; 
  vector<Type> psi = invlogit(logit_psi); 
  vector<Type> logit_p = X_p * beta_p; 
  if (S_p_n(0) > 0) logit_p += Z_p * z_p; 
  vector<Type> p = invlogit(logit_p); 
  
  // LIKELIHOOD 
  int i_psi = 0; 
  int i_y = 0; 
  for (int s = 0; s < nsites; ++s) {
    for (int k = 0; k < nocc(s); ++k) {
      if (totsite(i_psi) > 0) {
        nll -= log(psi(i_psi));
        for (int v = 0; v < nvisit(i_psi); ++v) {
          nll -= dbinom(y(i_y), Type(1.0), p(i_y), true);
          ++i_y;
        }
      } else {
        Type addllk = 0;
        for (int v = 0; v < nvisit(i_psi); ++v) {
          addllk += log(1 - p(i_y));
          ++i_y;
        }
        nll -= log(psi(i_psi) * exp(addllk) + 1 - psi(i_psi));
      }
      ++i_psi; 
    }
  }

  return nll;
}
