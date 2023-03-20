//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
functions {
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      real S = y[1];
      real I = y[2];
      real R = y[3];
      real N = x_i[1];
      
      real beta = theta[1];
      real gamma = theta[2];
      
      real dS_dt = -beta * I * S / N;
      real dI_dt =  beta * I * S / N - gamma * I;
      real dR_dt =  gamma * I;
      
      return {dS_dt, dI_dt, dR_dt};
  }
}
data {
  int M;                          //total number of data points
  int<lower=1> n_seasons;
  int<lower=1> n_weeks;
  int<lower=1> seasons[n_seasons];
  int<lower=1> seg_ind_start[n_seasons];
  int<lower=1> seg_ind_length[n_seasons];
  int<lower=1> seg_ind_max[n_seasons];
  real S0;
  real t0;
  real ts[M];                     //array of weeks
  int N;
  real ili[M];
  //int pos;
}
transformed data {
  real x_r[0];
  int x_i[1] = { N };
}
parameters {
  real<lower=0,upper=S0> rho;
  real<lower=0> beta;
  real<lower=0,upper=.1> I0;
  real<lower=0> phi_inv;
}
transformed parameters{
  real y0[3];
  real y[M, 3];
  real phi = 1. / phi_inv;
  real<lower=0> R0;
  real<lower=0> gamma;
  
  
  R0 = 1 - S0 - I0;
  gamma = rho*beta;
  //pos = 1;
  {
    real theta[2];
    theta[1] = beta;
    theta[2] = gamma;
    y0 = {S0, I0, R0};
    for (i in 1:n_seasons) {
      
      //print(dims(y[seg_ind_start[i]:seg_ind_max[i],]))
      y[seg_ind_start[i]:seg_ind_max[i],] =
        integrate_ode_rk45(sir, y0, t0, 
                            segment(ts, seg_ind_start[i], seg_ind_length[i]),
                            theta, x_r, x_i);
          //pos = pos + seg_ind_length[i];
    }
  }
}

model {
  //priors
  rho ~ normal(0.68, 0.08);
  beta ~ normal(.8, .3);
  I0 ~ normal(.005, .03);
  phi_inv ~ exponential(5);
  
  //sampling distribution
  //col(matrix x, int n) - The n-th column of matrix x. Here the number of infected people 
  for (i in 1:n_seasons) {
                        segment(ili, seg_ind_start[i], seg_ind_length[i])
                        ~ beta_proportion(col(
                            to_matrix(segment(y, 
                            seg_ind_start[i], seg_ind_length[i])), 2), 
                            phi);
  }
}
generated quantities {
  real r0 = beta / gamma;
  real recovery_time = 1 / gamma;
  //real yp[43];
  real pred_ili[2];
  real tp[2];
  real theta[2];
  real ypm[2,3];
  //real yp[2];
  //yp = y[1:43,2];
  
  tp = {22,23};
  
  
  theta[1] = beta;
  theta[2] = beta*rho;
  
  ypm = integrate_ode_rk45(sir, y0, t0, tp, theta, x_r, x_i);
  pred_ili = beta_proportion_rng(col(to_matrix(ypm) ,2), phi);
  //pred_ili = beta_proportion_rng(yp, phi);
  //for (i in 1:43) {
    //print(yp[i]);
  //  pred_ili[i] = beta_proportion_rng(yp[i], phi);
  //}
  
  //pred_ili = beta_proportion_rng(integrate_ode_rk45(sir, y0, t0, 
  //                          ts,
  //                          theta, x_r, x_i));
  
  
  
  
  
  
  
  
}