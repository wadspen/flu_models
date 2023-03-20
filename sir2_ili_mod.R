sir2_ili_mod = '
functions {
  
  real[] SIR(real t,
             real[] y,
             real[] theta,
             data real[] x_r,
             data int [] x_i) {
             
             real S = y[1];
             real I = y[2];
             real R = y[3];
             real N = x_i[1];
             
             real beta = theta[1];
             real gamma = theta[2];
             
             real dSdt = -beta * S * I / N;
             real dIdt = beta * S * I / N - gamma * I;
             real dRdt = gamma * I;
             
             return {dSdt, dIdt, dRdt};
             
         }
}

data {
int<lower=0> N;                                 // number of points
int<lower=0> n_seasons;                         // number of seasons
int<lower=0> seasons[N];                        // seasons
real<lower=0> ili[N];                           // all ili data
int<lower=0> x[N];                              // week
real<lower=0> n_weeks[43];       


int<lower=0> ps;                                // season to predict
real<lower=0,upper=1> S0;                       // initial susceptible proportion
real t0;                                        // week 0
}

transformed data {
  real x_r[0];
  int x_i[1] = { 1 };
}

parameters {


vector[n_seasons] beta_s;
vector[n_seasons] rho_s;
real<lower=0,upper=.1> I0_s[n_seasons];
real<lower=0> epsilon;


}



transformed parameters {

//matrix[n_seasons, n_compartments] y0;
real<lower=0,upper=1> R0;
vector[n_seasons] gamma_s;
//matrix[n_seasons, 2] theta;

for (i in 1:n_seasons) {
                gamma_s[i] = beta_s[i] * rho_s[i];
                R0 = 1 - S0 - I0_s[i];
                //y0[i,] = {S0, I0_s[i], R0};
                //theta[i,] = [beta_s[i], gamma_s[i]];
}



}

model {

for (i in 1:n_seasons) {
  beta_s[i] ~ gamma(2,.02);
  rho_s[i] ~ gamma(2,.02);
  I0_s[i] ~ normal(.05,.02);
}

epsilon ~ exponential(100);

//for (i in 1:N) ili[i] ~ normal(RK4SIR(x[i],beta_s[seasons[i]],
//                          rho_s[seasons[i]],.9,I0_s[seasons[i]],
//                          1,1)[x[i],2],
//                          epsilon);

for (i in 1:N) {
                 print("loop iteration: ", i);
                 ili[i] ~ normal(col(to_matrix(integrate_ode_rk45(SIR, 
                 {S0, I0_s[seasons[i]], 1 - S0 - I0_s[seasons[i]]}, 
                              t0, n_weeks, 
                              {beta_s[seasons[i]], gamma_s[seasons[i]]},
                              x_r, x_i)),2),
                              epsilon);
}
                          
                          
}

//generated quantities {
    //real ili_pred[42];
    //for (i in 1:42) {
      //real ili_pred[i] = normal_rng(integrate_ode_rk45(SIR, 
          //                    {S0, I0_s[ps], 1 - S0 - I0_s[ps]}, 
        //                      t0, {n_weeks[N]}, {beta_s[ps], gamma_s[ps]},
      //                        x_r, x_i)[1],
    //                          epsilon);
  //  }
//}
'