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
  vector[q] omega;
  vector[p] mu;
  real<lower=0> nu;
  real<lower=0> kappa;
  real<lower=0> lambda;
  real<lower=0> eta;
  real hs_1[p];
  real hs_2;
}

transformed parameters{
  vector[p] gamma;

  for (i in 1:p){
    gamma[i] = hs_1[i]*hs_2;
  }
}

model {

  tau ~ gamma(eta, kappa);
  sigma ~ gamma(lambda, nu);

  hs_1 ~ cauchy(0,1);
  hs_2 ~ cauchy(0,1);


  for (n in 1:p) {
    beta[n] ~ normal(gamma[n]*mu,gamma[n]*sigma*tau^(-1)+(1-gamma[n])*0.001);
  }

  y ~ normal(x*beta, sigma*tau^(-1));
}
