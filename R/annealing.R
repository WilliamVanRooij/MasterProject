
library(parallel)
library(LaplacesDemon)

require(MASS)
require(echoseq)
require(locus)
require(R2OpenBUGS)
require(coda)
require(lattice)
require(nor1mix)
require(plotly)

a <- 1

rm(list= ls())

RNGkind("L'Ecuyer-CMRG")
set.seed(123);

# ================================================================================================================================================ #
#                                                                                                                                                  #
#                                                                     TO DO                                                                        #
#                                                                                                                                                  #
# ================================================================================================================================================ #

# Check the Windows part of the function 

# ================================================================================================================================================ #
#                                                                                                                                                  #
#                                                                   VARIABLES                                                                      #
#                                                                                                                                                  #
# ================================================================================================================================================ #

n <- 300; 
p <- 500; p0 <- 5; 
d <- 1; d0 <- 1


iter <- 1


for (k in sample(1:1e3,iter)){
  
  set.seed(k)  
  
  seed <-  k;
  
  cor_type <- "autocorrelated"; 
  
  vec_rho <- runif(floor(p/10), min = 0.98, max = 0.99)
  
  vec_maf <- runif(p, 0.25,0.5)

  t_scheme <- 2
  
  t_final <- 5
  
  t_size <- 100
  
  nb_cpus <- 4;
  
  ind_d0 <-  sample(1:d, d0)
  
  ind_p0 <- c(3,13,17,23,43)
  
  p0_av <- 15
  
  user_seed <- sample(1:1e3, 100)
  
  vec_prob_sh <-  0.1 # proba that each SNP will be associated with another active phenotype
  
  max_tot_pve <-  0.4 # max proportion of phenotypic variance explained by the active SNPs
  
  list_snps <- generate_snps(n, p, cor_type, vec_rho, n_cpus = nb_cpus,
                             user_seed = seed, vec_maf = vec_maf)
  
  list_phenos <- generate_phenos(n, d,  user_seed = seed)
  
  dat_g <- generate_dependence(list_snps, list_phenos, ind_d0 = ind_d0,
                               ind_p0 = ind_p0, vec_prob_sh = vec_prob_sh,
                               family = "gaussian", max_tot_pve = max_tot_pve,
                               block_phenos = TRUE, user_seed = seed)
  
  
  mlocus <- function(fseed) {
    
    tau_vb_init <- 1 / apply(dat_g$phenos, 2, var)
    
    sig2_beta_vb_init <- 1 / rgamma(d, shape = 2, rate = 1 / tau_vb_init)
    
    gam_vb_init <-  matrix(rbeta(p * d, shape1 = 1, shape2 = d * (p - p0_av)/p0_av), nrow = p)
    
    
    mu_beta_vb_init <- matrix(rnorm(n = p * d , mean=0, sd = 1), nrow = p)
    
    list_init0 <-  set_init(d,p, gam_vb = gam_vb_init, mu_beta_vb = mu_beta_vb_init, 
                            sig2_beta_vb = sig2_beta_vb_init, tau_vb = tau_vb_init)
    
    vb_g <- locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0_av, link = "identity", user_seed = fseed, list_init = list_init0, full_output = TRUE, anneal = NULL)
    
    return(list(locus = vb_g, gam_init = gam_vb_init,mu = mu_beta_vb_init, sig = sig2_beta_vb_init))
  }
  
  mlocus_anneal <- function(fseed) {
    
    tau_vb_init <- 1 / apply(dat_g$phenos, 2, var)
    
    sig2_beta_vb_init <- 1 / rgamma(d, shape = 2, rate = 1 / tau_vb_init)
    
    gam_vb_init <-  matrix(rbeta(p * d, shape1 = 1, shape2 = d * (p - p0_av)/p0_av), nrow = p)
    
    
    mu_beta_vb_init <- matrix(rnorm(n = p * d , mean=0, sd = 1), nrow = p)
    
    list_init0 <-  set_init(d,p, gam_vb = gam_vb_init, mu_beta_vb = mu_beta_vb_init, 
                            sig2_beta_vb = sig2_beta_vb_init, tau_vb = tau_vb_init)
    
    vb_g <- locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0_av, link = "identity", user_seed = fseed, list_init = list_init0, full_output = TRUE, anneal = c(t_scheme,t_final,t_size))

    return(list(locus = vb_g, gam_init = gam_vb_init,mu = mu_beta_vb_init, sig = sig2_beta_vb_init))
  }
  
  
  # MAC
  m_vb_g <- mclapply(user_seed, mlocus, mc.cores = nb_cpus)
  m_vb_g_a <- mclapply(user_seed, mlocus_anneal, mc.cores = 1)
  
  if(TRUE) {
    
    
    out <- 0
    lb_exp <- 0
    elbo <- NULL
    out_a <- 0
    lb_exp_a <- 0
    elbo_a <- NULL
    
    
    if(TRUE){
      for(i in c(1:length(user_seed))) {
        out <- out + m_vb_g[[i]]$locus$gam_vb*exp(m_vb_g[[i]]$locus$lb_opt)
        lb_exp <- lb_exp + exp(m_vb_g[[i]]$locus$lb_opt)
        elbo <- c(elbo, exp(m_vb_g[[i]]$locus$lb_opt))
      }
 
      for(i in c(1:length(user_seed))) {
        out_a <- out_a + m_vb_g_a[[i]]$locus$gam_vb*exp(m_vb_g_a[[i]]$locus$lb_opt)
        lb_exp_a <- lb_exp_a + exp(m_vb_g_a[[i]]$locus$lb_opt)
        elbo_a <- c(elbo_a, exp(m_vb_g_a[[i]]$locus$lb_opt))
      }
    }
    
    sum_elbo <- get_p_m_y(elbo)
    out <- out / sum(sum_elb0)
    out_a <- out_a / lb_exp_a
    
  }
  
  single_vb_g <-locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0_av, link = "identity", user_seed = seed, verbose = FALSE, save_hyper=TRUE, save_init = TRUE,full_output=TRUE, anneal=NULL)
  single_vb_g_a <-locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0_av, link = "identity", user_seed = seed, verbose = FALSE, save_hyper=TRUE, save_init = TRUE,full_output=TRUE, anneal=c(t_scheme,t_final,t_size))
  
  if(FALSE) {
    plot(out_a[1:50], main='Probabilities of link between a phenotype and SNPs, comparison between m_locus and annealing',type='h',lwd=3,lend=1, ylim = c(0,1))
    points(single_vb_g$gam_vb[1:50],lwd=2,lend=1,type="h",col="red")
    points(ind_p0, out[ind_p0], col='black')
    points(ind_p0, single_vb_g$gam_vb[ind_p0], col = "red")
  }
  
  
  if(TRUE){
    par(mfrow=c(2,2))
    
    plot(single_vb_g$gam_vb[1:50],type='h',lwd=8,lend=1,xlab='',col='#a4a4a4', xaxt='n',main ="Probability of association - single LOCUS", ylab='')
    points(ind_p0, single_vb_g$gam_vb[ind_p0],col='red', type='h', lend=1,lwd=8)
    
    plot(single_vb_g_a$gam_vb[1:50],type='h',lwd=8,lend=1,xlab='',col='#a4a4a4', xaxt='n',main ="Probability of association - Annealed single LOCUS", ylab='')
    points(ind_p0, single_vb_g_a$gam_vb[ind_p0],col='red', type='h', lend=1,lwd=8)
    
    plot(out[1:50],type='h',lwd=8,lend=1,xlab='',col='#a4a4a4', xaxt='n',main ="Probability of association - Multiple LOCUS", ylab='')
    points(ind_p0, out[ind_p0], col = "red",type='h',lend=1,lwd=8)

    plot(out_a[1:50],type='h',lwd=8,lend=1,xlab='',col='#a4a4a4', xaxt='n',main ="Probability of association - Annealed multiple LOCUS", ylab='')
    points(ind_p0, out_a[ind_p0], col = "red",type='h',lend=1,lwd=8)
    
    par(mfrow=c(1,1))
  }
  
  
}




