sir_ili_mod = "
functions {

   real dSdt(real t, real S, real I, real beta, int N) {
     return -beta*S*I/N;       
   }
   
   real dIdt(real t, real S, real I, real beta, real gamma, int N) {
     return beta*S*I/N - gamma*I;
   }
   
   real dRdt(real t, real I, real gamma) {
     return gamma*I;
   }
   
   matrix RK4SIR(int n, real beta, real rho, real S0, real I0, real dt,
                    int N) {
                    
     real gamma;
     real R0;
     vector[n + 1] S;
     vector[n + 1] I;
     vector[n + 1] R;
     real Si;
     real Ii;
     //real Ri;
     real S_k1;
     real S_k2;
     real S_k3;
     real S_k4;
     real I_k1;
     real I_k2;
     real I_k3;
     real I_k4;
     //real R_k1;
     //real R_k2;
     //real R_k3;
     //real R_k4;
     matrix[n + 1,3] SIR;
     
     gamma = beta*rho;
     R0 = N - S0 - I0;
     
     S[1] = S0;
     I[1] = I0;
     R[1] = R0;
     
     for (i in 1:n) {
        Si = S[i];
        Ii = I[i];
        // R[i];
        
        S_k1 = dSdt(i, Si, Ii, beta, N);
        I_k1 = dIdt(i, Si, Ii, beta, gamma, N);
        //     R_k1 = dRdt(i, Ii);
    
        S_k2 = dSdt(i + dt / 2, Si + dt / 2 * S_k1, Ii + dt / 2 * I_k1, beta, 
                    N);
        I_k2 = dIdt(i + dt / 2, Si + dt / 2 * S_k1, Ii + dt / 2 * I_k1,
                 beta, gamma, N);
        //     R_k2 = dRdt(i + dt / 2, Ii + dt / 2 * I_k1);
    
        S_k3 = dSdt(i + dt / 2, Si + dt / 2 * S_k2, Ii + dt / 2 * I_k2, beta, 
                    N);
        I_k3 = dIdt(i + dt / 2, Si + dt / 2 * S_k2, Ii + dt / 2 * I_k2, 
                 beta, gamma, N);
        //     R_k3 = dRdt(i + dt / 2, Ii + dt / 2 * I_k2);
    
        S_k4 = dSdt(i + dt, Si + dt * S_k3, Ii + dt * I_k3, beta, N);
        I_k4 = dIdt(i + dt, Si + dt * S_k3, Ii + dt * I_k3, 
                 beta, gamma, N);
        //     R_k4 = dRdt(i + dt, Ii + dt * I_k3);
    
        S[i + 1] = Si + dt / 6 * (S_k1 + 2 * S_k2 + 2 * S_k3 + S_k4);
        I[i + 1] = Ii + dt / 6 * (I_k1 + 2 * I_k2 + 2 * I_k3 + I_k4);
        //     R[i + 1] = Ri + dt / 6 * (R_k1 + 2 * R_k2 + 2 * R_k3 + R_k4);
     }
     
     R = N - S - I;
     SIR[,1] = S;
     SIR[,2] = I;
     SIR[,3] = R;
     
     return SIR;
   }
}

data {
int<lower=0> N;                                 // number of points
int<lower=0> n_seasons;                         // number of seasons
int<lower=0> seasons[N];                        // seasons
real<lower=0> ili[N];                           // all ili data
int<lower=0> x[N];                              // week


int<lower=0> ps;                                // season to predict
int<lower=0> last_week;
}
parameters {


vector[n_seasons] beta_s;
vector[n_seasons] rho_s;
real<lower=0,upper=.1> I0_s[n_seasons];
real<lower=0> epsilon;


}

transformed parameters {


}

model {

for (i in 1:n_seasons) {
  beta_s[i] ~ gamma(2,.02);
  rho_s[i] ~ gamma(2,.02);
  I0_s[i] ~ normal(.05,.02);
}

epsilon ~ exponential(1);

for (i in 1:N) ili[i] ~ normal(RK4SIR(x[i],beta_s[seasons[i]],
                          rho_s[seasons[i]],.9,I0_s[seasons[i]],
                          1,1)[x[i],2],
                          epsilon);
                          
                          
}

generated quantities {
    real ili_pred[22];
    for (i in 1:22) {
      ili_pred[i] = normal_rng(RK4SIR(i,beta_s[ps],
                          rho_s[ps],.9,I0_s[ps],
                          1,1)[i,2]
                          ,epsilon);
    }
}
"