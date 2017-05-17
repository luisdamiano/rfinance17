data {
  int<lower=0> N;       # Number of assets
  int<lower=0> K;       # Number of states
  int<lower=0> T;       # Length of series
  real y[N, T];         # Series
}

parameters {
  // GARCH
  positive_ordered[K] sigma[N]; # Ordering prevents in-chain switching
  real<lower=0.00000001> nu;
  
  real<lower=0.00000001, upper=1> alpha1[N, K];
  real<lower=0, upper=1-alpha1[1, 1]> beta11;
  real<lower=0, upper=1-alpha1[2, 1]> beta21;
  real<lower=0, upper=1-alpha1[3, 1]> beta31;
  real<lower=0, upper=1-alpha1[4, 1]> beta41;
  real<lower=0, upper=1-alpha1[5, 1]> beta51;
  real<lower=0, upper=1-alpha1[1, 2]> beta12;
  real<lower=0, upper=1-alpha1[2, 2]> beta22;
  real<lower=0, upper=1-alpha1[3, 2]> beta32;
  real<lower=0, upper=1-alpha1[4, 2]> beta42;
  real<lower=0, upper=1-alpha1[5, 2]> beta52;

  // Regime switching
  real<lower=0, upper=1> p0[N] ;     // initial prob for K = 1
  real<lower=0, upper=1> TP[N, K] ;  // transition probs of staying in group K (tau)
}

transformed parameters {
  // GARCH
  real<lower=0.00000001> sigma_t[N, K, T];
  real<lower=0.00000001> alpha0[N, K];
  real<lower=0, upper=1> beta[N, K];
  vector<lower=0, upper=1>[K] pi_t[N, T];

  // Regime switching
  real F[N, T]; //filter forwards group membership prob
  real p1;
  real p2;

  // GARCH
  beta[1, 1] = beta11;
  beta[2, 1] = beta21;
  beta[3, 1] = beta31;
  beta[4, 1] = beta41;
  beta[5, 1] = beta51;
  beta[1, 2] = beta12;
  beta[2, 2] = beta22;
  beta[3, 2] = beta32;
  beta[4, 2] = beta42;
  beta[5, 2] = beta52;

  for (n in 1:N) {      # For each asset n = 1, ..., N
    for (k in 1:K) {    # For each state k = 1, ..., K
      sigma_t[n, k, 1] = sigma[n][k]; # Sets the cond vol for the first period to the uncond vol
      # alpha0[n, k] = pow(sigma[n][k], 2) * (1 - alpha1[n, k] - beta[n, k]);
      alpha0[n, k] = exp(2 * log(sigma[n][k]) + log(1 - alpha1[n, k] - beta[n, k]));
      
      for (t in 2:T) {  # For each period t = 1, ..., T
        sigma_t[n, k, t] = sqrt(alpha0[n, k]
                           + alpha1[n, k] * pow(y[n, t-1], 2)
                           + beta[n, k] * pow(sigma_t[n, k, t-1], 2));
      }
    }
  }
  
  // FORWARD ALGORITHM
  for (n in 1:N) {
    F[n, 1] = p0[n];
    pi_t[n, 1][1] = F[n, 1];
    pi_t[n, 1][2] = 1 - F[n, 1];
  }
  
  for (t in 1:T){
    for (n in 1:N) {
      //1 update prior using data (likelihood) exp
      p1 = exp(normal_lpdf(y[n, t] | 0, sigma_t[n, 1, t])) * F[n, t];
      p2 = exp(normal_lpdf(y[n, t] | 0, sigma_t[n, 2, t])) * (1-F[n, t]);
      F[n, t] = p1/(p1+p2);

      //2 probability - forward filtered one step ahead
      if (t != T) {
        p1 = F[n, t]*TP[n, 1] + (1-F[n, t])*(1-TP[n, 2]);
        p2 = F[n, t]*(1-TP[n, 1]) + (1-F[n, t])*TP[n, 2];
        F[n, t+1] = p1/(p1+p2);
        pi_t[n, t+1][1] = F[n, t+1];
        pi_t[n, t+1][2] = 1 - F[n, t+1];
      }
    }
  }
}

model {
  # Priors
  nu ~ normal(0,  5);
  for(n in 1:N) {
    sigma[n] ~ normal(0,  nu);      # sigma prior needed to easy stuck chains 
                                    # see https://groups.google.com/forum/#!topic/stan-users/xEG18UzoCGo
  }

  # Likelihood
  for(n in 1:N) {
    for(k in 1:K) {
      real tmp_accum[T];
      for(t in 1:T) {
        tmp_accum[t] = pi_t[n, t][k] * normal_lpdf(y[n, t] | 0, sigma_t[n, k, t]);
      }
      
      target += log_sum_exp(tmp_accum);
    }
  }
}

generated quantities {

}
