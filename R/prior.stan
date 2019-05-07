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
  vector[p] gamma;
  vector[p] omega;
  vector[p] mu;
  real nu;
  real kappa;
  real lambda;
  real eta;
  real a;
  real b;

}
model {
  tau ~ gamma(eta, kappa);
  sigma ~ gamma(lambda, nu);
  omega ~ beta(a,b);
  gamma ~

  for (n in 1:p) {
    beta[n] ~ normal(gamma[n]*mu,gamma[n]*sigma*tau^(-1)+(1-gamma[n])*0.001);
  }

  y ~ normal(x*beta, sigma*tau^(-1));
}
