data {
  int<lower=1> N;         // number of observations
  int<lower=1> q;         // number of traits
  int<lower=1> p;         // number of predictors
  vector[N] y;            // observations, traits
  matrix[N,p] x;          // predictors, SNPs
}
parameters {
  real sigma;
  real tau;
  vector[p] beta;
  real omega;
  vector[p] mu;
  real<lower=0> nu;
  real<lower=0> kappa;
  real<lower=0> lambda;
  real<lower=0> eta;
  real<lower=0,upper=1> hs_1[p];
  real<lower=0,upper=1> hs_2;
}

transformed parameters{
  vector<lower=0,upper=1>[p] gamma;

  for (i in 1:p){
    gamma[i] = hs_1[i]*hs_2;
  }
}

model {

  tau ~ gamma(eta, kappa);
  sigma ~ gamma(lambda, nu);
  omega ~ beta(1,1);
  for (i in 1:p){
    hs_1[i] ~ normal(0,1);
  }
  hs_2 ~ normal(0,1);


  for (n in 1:p) {
    beta[n] ~ normal(gamma[n]*mu,gamma[n]*sigma*tau^(-1)+(1-gamma[n])*0.1);
  }

  y ~ normal(x*beta, sigma*tau^(-1));
}
