functions {
  vector log_2nd_poly_vector(vector c, vector x, vector y) {
    return(c[1] + c[2]*log(x).*log(x) + c[3]*log(y).*log(y) + c[4]*log(x) + c[5]*log(y) + c[6]*log(x).*log(y));
  }
  real log_2nd_poly_scalar(vector c, real x, real y) {
    return(c[1] + c[2]*log(x).*log(x) + c[3]*log(y).*log(y) + c[4]*log(x) + c[5]*log(y) + c[6]*log(x).*log(y));
  }
}

data {
  int Nt;
  int Ny;
  matrix[Nt, 3] EEUt;
  matrix[Ny, 3] EEUy;
  int Lt;
  int Ly;
  vector[Lt] Yt;
  vector[Ly] Yy;
  vector[6] Ct;
  vector[6] Cy;
  matrix[4, 2] Bounds;
}

parameters {
  real<lower=0, upper=0.1> sigma;
  real<lower=0, upper=0.1> lambda;
  real<lower=Bounds[1,1], upper=Bounds[1,2]> e1;
  real<lower=Bounds[2,1], upper=Bounds[2,2]> e2;
  real<lower=Bounds[3,1], upper=Bounds[3,2]> ut;
  real<lower=Bounds[4,1], upper=Bounds[4,2]> uy;
}

transformed parameters {
  vector[Nt] Ut;
  vector[Ny] Uy;
  vector<lower=min(EEUt[, 1]), upper=max(EEUt[, 1])>[Nt] E1t;
  vector<lower=min(EEUy[, 1]), upper=max(EEUy[, 1])>[Ny] E1y;
  vector<lower=min(EEUt[, 2]), upper=max(EEUt[, 2])>[Nt] E2t;
  vector<lower=min(EEUy[, 2]), upper=max(EEUy[, 2])>[Ny] E2y;
  
  E1t = EEUt[, 1];
  E1y = EEUy[, 1];
  E2t = EEUt[, 2];
  E2y = EEUy[, 2];
  Ut =  EEUt[, 3];
  Uy =  EEUy[, 3];
}

model {
  //  lambda ~ normal(0, 0.01);
  lambda ~  gamma(2, 0.1);
  sigma ~ normal(0, 0.01);
  
  Yt ~ normal(ut, lambda);
  Yy ~ normal(uy, lambda);
  // Yt ~ normal(ut, 1/sqrt(lambda));
  // Yy ~ normal(uy, 1/sqrt(lambda));
  
  Ut ~ normal(log_2nd_poly_vector(Ct, E1t, E2t), sigma);
  Uy ~ normal(log_2nd_poly_vector(Cy, E1y, E2y), sigma);

  ut ~ normal(log_2nd_poly_scalar(Ct, e1, e2), sigma);
  uy ~ normal(log_2nd_poly_scalar(Cy, e1, e2), sigma);
}

generated quantities {
  vector[Lt+Ly] log_lik;
  real y_pred_t;
  real y_pred_y;

  for (ind in 1:Lt) {
    log_lik[ind] = normal_lpdf(Yt[ind] | ut, lambda);
  }
  for (ind in 1:Ly) {
    log_lik[ind+Lt] = normal_lpdf(Yy[ind] | uy, lambda);
  }

  y_pred_t = normal_rng(ut, lambda);
  y_pred_y = normal_rng(uy, lambda);
}
