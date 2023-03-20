asg_ili_2_disc_mod = "
functions {

   real lognormal_2_lpdf(real x, real mu, real sigma) {
     return 1/(x*sigma*sqrt(2*pi()))*exp((-(log(x) - mu)^2)/(2*sigma^2));       
   }

   real asg(row_vector theta, real x) {
    real beta1; 
    real beta2;
    real eta;
    real mu; 
    real sig1;
    real sig2;
    real ASG;
    
    beta1 = exp(theta[1]);
    beta2 = exp(theta[2]);
    eta = theta[3];
    mu = theta[4];
    sig1 = exp(theta[5]);
    sig2 = exp(theta[6]);
    
    ASG = ((beta1 + (eta-beta1)*exp(-((x-mu)^2)/(2*sig1^2)))*(x < mu) +
       (beta2 + (eta-beta2)*exp(-((x-mu)^2)/(2*sig2^2)))*(mu <= x));
    return ASG;
  }
}

data {
int<lower=0> N;                                 // number of points
int<lower=0> n_seasons;                         // number of seasons
int<lower=0> seasons[N];                        // seasons
real<lower=0> ili[N];                           // all ili data
int<lower=0> x[N];                              // week
int<lower=0> n_params;                          // number of parameters ASG
vector[n_params] m0;                            // prior means
matrix[n_params,n_params] C0;                   // prior sds

real nu;
real<lower=0> c;
real<lower=0> d;
int<lower=0> m;                                 // weeks in a season

vector[n_seasons] es;
matrix[n_seasons,n_seasons] ss;
int<lower=0> ps;                                // season to predict
int<lower=0> last_week;                         // last observed week
}
parameters {


vector[n_params] theta;
matrix[n_seasons,n_params] theta_s;
real<lower=0> tau;
real<lower=0> taumu;
real<lower=0> epsilon;
real<lower=0> epsilons[n_seasons];


row_vector[n_params] zeta;
corr_matrix[n_params] omega;
real gamma[m];

matrix[n_seasons,m] delta;
real<lower=0> sig_delt[n_seasons];
real<lower=0,upper=1> alpha[n_seasons];
real<lower=0> sig_alph;
real<lower=0> a_delt;
real<lower=0> b_delt;

}

transformed parameters {
matrix[n_params,n_params] Sigma;


Sigma = quad_form_diag(omega,zeta) ;

}

model {

theta ~ multi_normal(m0,C0);
tau ~ gamma(2,2);
taumu ~ gamma(2,.2);
epsilon ~ exponential(1);
sig_alph ~ gamma(2,2);

zeta ~ student_t(4,c,d);
omega ~ lkj_corr(nu);

a_delt ~ gamma(5,1);
b_delt ~ gamma(1,10);
                    
for (i in 1:n_seasons) {
                          theta_s[i,] ~ multi_normal(theta,Sigma);
                          epsilons[i] ~ exponential(epsilon);
                          alpha[i] ~ normal(0,sig_alph^2);
                          sig_delt[i] ~ inv_gamma(a_delt,b_delt);
}

gamma[m] ~ normal(0,1/tau);
delta[,m] = theta_s[,2];

for (i in 1:(m-1)) {gamma[i] ~ normal(gamma[i+1],1/taumu);
                    delta[,i] ~ multi_normal(alpha*delta[,i+1],
                                             diag_matrix(1/sig_delt));
}

for (i in 1:N) ili[i] ~ ili[i] ~ beta_proportion(
                                    inv_logit(asg(theta_s[seasons[i],],x[i]) + 
                                    gamma[x[i]] + delta[seasons[i],x[i]]),
                                    epsilons[seasons[i]]);
                          

}

generated quantities {
    real ili_pred[4];
    for (i in 1:4) {
      //ili_pred[i] = normal_rng(asg(theta_s[ps,],x[last_week + i]) + 
      //      gamma[last_week + i],epsilons[ps]);
            
      ili_pred[i] = beta_proportion_rng(
                    inv_logit(asg(theta_s[ps,],last_week + i) + 
                                    gamma[last_week + i] + delta[ps,
                                        last_week + i]),
                                    epsilons[ps]          
      )
        
    }
}
"

# library(rstan)
mo <- stan_model(model_code = asg_ili_2_disc_mod)
