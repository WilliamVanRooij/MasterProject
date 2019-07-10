library(parallel)

require(echoseq)
require(locus)
require(ROCR)

rm(list= ls())

RNGkind("L'Ecuyer-CMRG")
seed <- 166
set.seed(seed);


n <- 300; # Number of observations
p <- 500; p0 <- 50; #Number of SNPs ; Number of active SNPs
d <- 1; d0 <- 1 # Number of traits ; Number of active traits

max_tot_pve <-  0.8 # max proportion of phenotypic variance explained by the active SNPs

min_rho <- 0.5; max_rho <- 0.99; # Minimum and maximum correlation between the SNPs 
iter <- 50

anneal <- c(1, 2, 10)

cor_type <- "autocorrelated"; 


# Definition of useful functions

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

make_ld_plot <- function(X, meas) {
  
  stopifnot(meas %in% c("r", "D'"))
  
  require(LDheatmap)
  require(chopsticks)
  
  colnames(X)<- paste(1:ncol(X), "  ", sep="")
  require(snpStats)
  gX <- as(X, "SnpMatrix")
  
  cat("LD plot display:\n")
  ld <- LDheatmap(gX, flip=TRUE, name="", title=NULL, LDmeasure = meas,
                  add.map= T, geneMapLocation = 0.01, geneMapLabelX=1000)
}


# Definiton of fucntion for Averaged LOCUS

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

# Definition of funciton for Averaged annealed LOCUS

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

  
  c_pred <- NULL
  c_lab <- NULL
  
  c_pred_m <- NULL
  c_lab_m <- NULL
  
  single_pred <-  NULL
  single_lab <-  NULL
  
  c_pred_a <- NULL
  c_lab_a <- NULL
  
  single_pred_a <-  NULL
  single_lab_a <-  NULL


for(seed in sample(1:1e3,iter)){
  
  set.seed(seed)

  nb_cpus <- 4;
  
  ind_d0 <-  sample(1:d, d0)
  
  ind_p0 <- sample(1:p, p0)
  
  p0_av <- 50
  
  vec_rho <- runif(floor(p/10), min = min_rho, max = max_rho) # Autocorrelations between the SNPs
  
  vec_maf <- runif(p, 0.4, 0.5)

  vec_prob_sh <-  0.05 # proba that each SNP will be associated with another active phenotype, Not used when number of traits is one.
  
  list_snps <- generate_snps(n, p, cor_type, vec_rho, n_cpus = nb_cpus,
                             user_seed = seed, vec_maf = vec_maf)
  
  list_phenos <- generate_phenos(n, d,  user_seed = seed)
  
  dat_g <- generate_dependence(list_snps, list_phenos, ind_d0 = ind_d0,
                               ind_p0 = ind_p0, vec_prob_sh = vec_prob_sh,
                               family = "gaussian", max_tot_pve = max_tot_pve,
                               block_phenos = TRUE, user_seed = seed)
  
  
  user_seed <- sample(1:1e3, 100)
  
    m_vb_g <- mclapply(user_seed, mlocus, mc.cores = nb_cpus)
    m_vb_g_a <- mclapply(user_seed, a_mlocus, mc.cores = nb_cpus) 
    
    elbo <- NULL
    gam <- NULL

    elbo_a <- NULL
    gam_a <- NULL

    out_m <- 0

    

      for(i in c(1:length(user_seed))) {
        gam <- rbind(gam, as.vector(m_vb_g[[i]]$gam_vb))
        elbo <- c(elbo, m_vb_g[[i]]$lb_opt)

        gam_a <- rbind(gam_a, as.vector(m_vb_g_a[[i]]$gam_vb))
        elbo_a <- c(elbo_a, m_vb_g_a[[i]]$lb_opt)
        
        out_m <- out_m + m_vb_g[[i]]$gam_vb
    }
    
    vec_w_part <- get_p_m_y(elbo)
    out <- colSums(sweep(gam, 1, vec_w_part, "*")) # Averaged LOCUS gammas
    
    vec_w_part_a <- get_p_m_y(elbo_a)
    out_a <- colSums(sweep(gam_a, 1, vec_w_part_a, "*")) # Averaged annealed LOCUS gammas
    
    single_vb_g <-locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0_av, link = "identity", user_seed= seed, verbose = FALSE, save_hyper=TRUE, anneal = NULL) # LOCUS gammas
    single_vb_g_a <-locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0_av, link = "identity", user_seed = seed, verbose = FALSE, save_hyper=TRUE, anneal = anneal) # Annealed LOCUS gammas

    out_m <-  out_m/length(user_seed) # Averaged LOCUS (Equal weights) gammas

    #Define the performances
    
    single_pred <- cbind(single_pred, single_vb_g$gam_vb)
    single_lab <- cbind(single_lab, c(1:500) %in% ind_p0)
    
    c_pred <- cbind(c_pred, out)
    c_lab <- cbind(c_lab, c(1:500) %in% ind_p0)
    
    single_pred_a <- cbind(single_pred_a, single_vb_g_a$gam_vb)
    single_lab_a <- cbind(single_lab_a, c(1:500) %in% ind_p0)
    
    c_pred_a <- cbind(c_pred_a, out_a)
    c_lab_a <- cbind(c_lab_a, c(1:500) %in% ind_p0)
    
    c_pred_m <- cbind(c_pred_m, out_m)
    c_lab_m <- cbind(c_lab_m, c(1:500) %in% ind_p0)
}

  pred_m_locus <- prediction(c_pred, c_lab)
  pred_s_locus <- prediction(single_pred, single_lab)

  perf_m_locus <- performance(pred_m_locus, "tpr","fpr")
  perf_s_locus <- performance(pred_s_locus, "tpr","fpr")
  
  pred_m_locus_a <- prediction(c_pred_a, c_lab_a)
  pred_s_locus_a <- prediction(single_pred_a, single_lab_a)
  
  perf_m_locus_a <- performance(pred_m_locus_a, "tpr","fpr")
  perf_s_locus_a <- performance(pred_s_locus_a, "tpr","fpr")
  
  pred_m_locus_m <- prediction(c_pred_m, c_lab_m)
  perf_m_locus_m <- performance(pred_m_locus_m, "tpr","fpr")
  
  # Plot ROC curves
  
  par(pty="s")
  plot(perf_m_locus,avg="vertical",spread.estimate="stderror",spread.scale=2,col='orange',lwd=2, main=expression(paste("ROC Curves comparison, ",p[0]," = 50, Max Tot. PVE = 0.8", sep="")),xlim=c(0,0.2))
  plot(perf_s_locus,avg="vertical",spread.estimate="stderror",spread.scale=2,col='blue', lwd=2, add=T)
  plot(perf_m_locus_a,avg="vertical",spread.estimate="stderror",spread.scale=2,col='red', lwd=2, add=T)
  plot(perf_s_locus_a,avg="vertical",spread.estimate="stderror",spread.scale=2,col='green', lwd=2, add=T)
  plot(perf_m_locus_m,avg="vertical",spread.estimate="stderror",spread.scale=2,col='purple',lwd=2,add=T)
  legend(0.05,0.25, c("LOCUS","Annealed LOCUS", "Averaged LOCUS","Averaged annealed LOCUS","Averaged LOCUS (Equal weights)"), col=c('blue', 'green','orange', 'red','purple'),lwd=2)

  par(pty="m")
  