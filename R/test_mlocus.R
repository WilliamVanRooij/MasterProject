library(parallel)
require(echoseq)
require(locus)

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

seed <- 123; 
set.seed(seed);

n <- 200; 
p <- 500; p0 <- 5; 
d <- 1; d0 <- 1

cor_type <- "autocorrelated"; 

vec_rho <- runif(floor(p/20), min = 0.95, max = 0.99)

nb_cpus <- 4;

ind_d0 <-  sample(1:d, d0)

ind_p0 <- sample(1:p, p0)

user_seed <- sample(1:1e3, 100)

vec_prob_sh <-  0.05 # proba that each SNP will be associated with another active phenotype

max_tot_pve <-  0.5 # max proportion of phenotypic variance explained by the active SNPs

gam_vb_init <-  matrix(rbeta(p * d, shape1 = 1, shape2 = 4*d-1), nrow = p)

mu_beta_vb_init <- matrix(rnorm(p * d), nrow = p)



list_snps <- generate_snps(n, p, cor_type, vec_rho, n_cpus = nb_cpus,
                           user_seed = seed)

list_phenos <- generate_phenos(n, d, user_seed = seed)

dat_g <- generate_dependence(list_snps, list_phenos, ind_d0 = ind_d0,
                             ind_p0 = ind_p0, vec_prob_sh = vec_prob_sh,
                             family = "gaussian", max_tot_pve = max_tot_pve, 
                             block_phenos = TRUE, user_seed = seed)



tau_vb_init <- 1 / apply(dat_g$phenos, 2, var)

sig2_beta_vb_init <- 1 / rgamma(d, shape = 2, rate = 1 / tau_vb_init)



list_init <-  set_init(d,p, gam_vb = gam_vb_init, mu_beta_vb = mu_beta_vb_init, 
                       sig2_beta_vb = sig2_beta_vb_init, tau_vb = tau_vb_init)

params <- list(Y = dat_g$phenos, X = dat_g$snps, p0_av = 100, link = "identity", list_init <- list_init);


# ================================================================================================================================================ #
#                                                                                                                                                  #
#                                                                  FUNCTIONS                                                                       #
#                                                                                                                                                  #
# ================================================================================================================================================ #


mlocus <- function(fseed) {
  vb_g <- locus(Y = params$Y, X=params$X, p0_av = params$p0_av, link = params$link, user_seed = fseed,list_init = params$list_init)
  return(vb_g$gam_vb) # if d>1 and ind_d0 sampled randomly, this response may be inactive (i.e. no SNP associated with it, so there will be nothing to see...)
}

mac = (Sys.info()['sysname'] != "Windows")

if(!mac) {
  
  cores <- detectCores()
  cl <- makeCluster(cores)
  out <- clusterApply(cl=cl, x=user_seed , fun=mlocus) # are you sure you still get a vector here? maybe see below?
  plot(out, main='Probabilities of link between a phenotype and SNPs',type='h',lwd=2,lend=1, ylim = c(0,1))
  points(ind_p0, out[ind_p0], col = "red")

  }

# MAC
if(mac) {
  
  m_vb_g <- mclapply(user_seed, mlocus, mc.cores = nb_cpus) 
  
  out <- Reduce('+',  m_vb_g) / length(user_seed) # a bit more compact and no "hard coded" numbers
  plot(out, main='Probabilities of link between a phenotype and SNPs',type='h',lwd=2,lend=1, ylim = c(0,1))
  points(ind_p0, out[ind_p0], col = "red")
  # sum <- 0*(1:p)
  
}


# id <- 1
# plot(m_vb_g[[id]], main='Probabilities of link between a phenotype and SNPs',type='h',lwd=2,lend=1, ylim = c(0,1))
# points(ind_p0, out[ind_p0], col = "red")
# Next steps:
# implement weights in average
# compare with results from a single seed
# assess the performance with ROC curves, e.g., using the ROCR package
# think of a 2D example where we can visualize the local modes, e.g., 
# with two highly correlated predictors. See if the multiple-seed algorithm 
# is about to explore all the local modes (inspiration in simulation of Rockova et al. papers).
