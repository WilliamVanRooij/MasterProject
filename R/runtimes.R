library(parallel)

require(echoseq)
require(locus)
require(ROCR)

rm(list= ls())

RNGkind("L'Ecuyer-CMRG") 
seed <- 166
set.seed(seed);


n <- 300; 
p <- 500; p0 <- 15; # Number of SNPs ; Number of active SNPs
d <- 1; d0 <- 1 # Number of traits ; Number of active traits 

iter <- 50 # Number of replications 

anneal <- c(1, 2, 10) # Determination of the annealing parameters

min_rho <- 0.98; max_rho <- 0.99 # Minimum and maximum correlation between SNPs

runtime_s_a <- runtime_s <- runtime_m_a <- runtime_m <- NULL

#Definition of useful functions

log_sum_exp_ <- function(x) { # avoid numerical underflow or overflow
  if ( max(abs(x)) > max(x) )
    offset <- min(x)
  else
    offset <- max(x)
  log(sum(exp(x - offset))) + offset
}

get_p_m_y <- function(vec_elbo) {
  
  exp(vec_elbo - log_sum_exp_(vec_elbo))
}


# Definition of the function for averaged LOCUS

mlocus <- function(fseed) {
  
  tau_vb_init <- 1 / apply(dat_g$phenos, 2, var)
  
  sig2_beta_vb_init <- 1 / rgamma(d, shape = 2, rate = 1 / tau_vb_init)
  
  gam_vb_init <-  matrix(rbeta(p * d, shape1 = 1, shape2 = 4*d-1), nrow = p)
  
  mu_beta_vb_init <- matrix(rnorm(n = p * d , mean=0, sd = 1), nrow = p)
  
  list_init0 <-  set_init(d,p, gam_vb = gam_vb_init, mu_beta_vb = mu_beta_vb_init, 
                          sig2_beta_vb = sig2_beta_vb_init, tau_vb = tau_vb_init)
  
  vb_g <- locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0_av, link = "identity", user_seed = fseed, list_init = list_init0, save_hyper=TRUE, anneal = NULL)
  return(vb_g)
}


# Definition of the function for averaged annealed LOCUS

a_mlocus <- function(fseed) {
  
  tau_vb_init <- 1 / apply(dat_g$phenos, 2, var)
  
  sig2_beta_vb_init <- 1 / rgamma(d, shape = 2, rate = 1 / tau_vb_init)
  
  gam_vb_init <-  matrix(rbeta(p * d, shape1 = 1, shape2 = 4*d-1), nrow = p)
  
  mu_beta_vb_init <- matrix(rnorm(n = p * d , mean=0, sd = 1), nrow = p)
  
  list_init0 <-  set_init(d,p, gam_vb = gam_vb_init, mu_beta_vb = mu_beta_vb_init, 
                          sig2_beta_vb = sig2_beta_vb_init, tau_vb = tau_vb_init)
  
  vb_g <- locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0_av, link = "identity", user_seed = fseed, list_init = list_init0, save_hyper=TRUE, anneal = anneal)
  return(vb_g)
}

for(seed in sample(1:1e3,iter)){
  
  set.seed(seed)
  
  cor_type <- "autocorrelated"; 
  
  vec_rho <- runif(floor(p/10), min = min_rho, max = max_rho)

  nb_cpus <- 4;
  
  ind_d0 <-  sample(1:d, d0)
  
  ind_p0 <- sample(1:p, p0)
  
  p0_av <- 1
  
  vec_maf <- runif(p, 0.4, 0.5)

  vec_prob_sh <-  0.05 # proba that each SNP will be associated with another active phenotype
  
  max_tot_pve <-  0.5 # max proportion of phenotypic variance explained by the active SNPs
  
  list_snps <- generate_snps(n, p, cor_type, vec_rho, n_cpus = nb_cpus,
                             user_seed = seed, vec_maf = vec_maf)
  
  list_phenos <- generate_phenos(n, d,  user_seed = seed)
  
  dat_g <- generate_dependence(list_snps, list_phenos, ind_d0 = ind_d0,
                               ind_p0 = ind_p0, vec_prob_sh = vec_prob_sh,
                               family = "gaussian", max_tot_pve = max_tot_pve,
                               block_phenos = TRUE, user_seed = seed)
  
  
  user_seed <- sample(1:1e3, 100) # Determination of the seeds for the random initalisations
  
  
  time0_m <- proc.time()
  m_vb_g <- mclapply(user_seed, mlocus, mc.cores = nb_cpus)
  
  elbo <- NULL
  gam <- NULL

  for(i in c(1:length(user_seed))) {
    gam <- rbind(gam, as.vector(m_vb_g[[i]]$gam_vb))
    elbo <- c(elbo, m_vb_g[[i]]$lb_opt)
  }
  
  vec_w_part <- get_p_m_y(elbo)
  out <- colSums(sweep(gam, 1, vec_w_part, "*"))
  runtime_m <- c(runtime_m, as.numeric(proc.time() - time0_m)[3]);
  
  
  time0_m_a <- proc.time()
  m_vb_g_a <- mclapply(user_seed, a_mlocus, mc.cores = nb_cpus)
  
  elbo_a <- NULL
  gam_a <- NULL

  for(i in c(1:length(user_seed))) {
    gam_a <- rbind(gam_a, as.vector(m_vb_g_a[[i]]$gam_vb))
    elbo_a <- c(elbo_a, m_vb_g_a[[i]]$lb_opt)
  }
  
  vec_w_part_a <- get_p_m_y(elbo_a)
  out_a <- colSums(sweep(gam_a, 1, vec_w_part_a, "*"))
  runtime_m_a <- c(runtime_m_a, as.numeric(proc.time() - time0_m_a)[3])
  
  
  time0_s_a <- proc.time()
  single_vb_g_a <-locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0_av, link = "identity", user_seed = seed, verbose = FALSE, save_hyper=TRUE, anneal = anneal)
  runtime_s_a <- c(runtime_s_a, as.numeric(proc.time() - time0_s_a)[3])
  
  time0_s <- proc.time()
  single_vb_g <-locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0_av, link = "identity", user_seed= seed, verbose = FALSE, save_hyper=TRUE, anneal = NULL)
  runtime_s <- c(runtime_s, as.numeric(proc.time() - time0_s)[3])
  

}

par(mfrow=c(1,1))

boxplot(runtime_s, runtime_s_a, runtime_m, runtime_m_a,col=c("blue","green","orange","red"),lend=1, main="Running times of the four methods (in seconds)", xaxt='n',xlab="",ylab="Runtimes")
axis(1,at=1:4,labels=c("LOCUS","Annealed LOCUS", "Averaged LOCUS","Averaged annealed LOCUS"))

