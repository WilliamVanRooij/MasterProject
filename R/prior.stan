data {
  int<lower=1> N;         // number of observations
  int<lower=1> q;         // number of traits
  int<lower=1> p;         // number of predictors
  matrix[N,1] y;          // observations, traits
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

}
model {
  tau ~ inv_gamma(eta, kappa);
  sigma ~ gamma(lambda, nu);
  for (n in 1:q) {
    beta[n] ~ normal(gamma[n]*mu,gamma[n]*sigma*tau+(1-gamma[n])*0.001);
  }
  y ~ multi_normal(beta'*x, sigma*tau);
}
