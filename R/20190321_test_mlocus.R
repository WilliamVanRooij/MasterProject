library(parallel)
require(echoseq)
require(locus)


i <- 1;

n <- 200; p <- 100; p0 <- 10; d <- 100; d0 <- 80 # just one response for now (d <- 1, d0 <- 1)

gam_vb <- 0*c(1:p0); # initalization not needed

cor_type <- "autocorrelated"; 

vec_rho <- runif(100, min = 0.25, max = 0.95)

nb_cpus <- 4; # Unfortunately, the number of cores available on Windows is only 1. But when computing on another computer, one can change the number of cpus used.

seed <- 123; set.seed(seed);

list_snps <- generate_snps(n, p, cor_type, vec_rho, n_cpus = nb_cpus,
                           user_seed = seed)

list_phenos <- generate_phenos(n, d, cor_type = cor_type, vec_rho = vec_rho,
                               n_cpus = nb_cpus, user_seed = seed)

dat_g <- generate_dependence(list_snps, list_phenos, ind_d0 = sample(1:d, d0),
                             ind_p0 = sample(1:p, p0), vec_prob_sh = 0.95, # too high (but if d = 1 not used)
                             family = "gaussian", max_tot_pve = 0.5, block_phenos = TRUE, user_seed = seed)

# vb_g <- locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0, link = "identity", user_seed = seed)



user_seed <- cbind(123,234,345,456,567,678,789,890,901,012) # c()

params <- list(Y = dat_g$phenos, X = dat_g$snps, p0_av = p0, link = "identity");

mlocus <- function(fseed) {
  vb_g <- locus(Y = params$Y, X=params$X, p0_av = params$p0_av, link = params$link, user_seed = fseed)
  return(vb_g$gam_vb[,1]) # if d>1 and ind_d0 sampled randomly, this response may be inactive (i.e. no SNP associated with it, so there will be nothing to see...)
}

mac = TRUE

## WINDOWS Sys.info['sysname']

if(!mac) {
  cores <- detectCores()

cl <- makeCluster(cores)

mean(unlist(clusterApply(cl=cl,x=user_seed,fun=mlocus))) # are you sure you still get a vector here? maybe see below?

}

# MAC
if(mac) {
  cores <- detectCores()
  
  m_vb_g <- mclapply(user_seed, mlocus, mc.cores = 1) # = n_cpus instead of 1? or = cores if using detectCores?
  
  out <- Reduce('+',  m_vb_g) / length(user_seed) # a bit more compact and no "hard coded" numbers
  plot(out) # what do you conclude?
  # sum <- 0*(1:p)
  
  # for(i in (1:10)){
  #    sum <- sum + m_vb_g[[i]]
  # }
  # plot(sum/10)
  
}

# Next steps:
# compare with results from a single seed
# assess the performance with ROC curves, e.g., using the ROCR package
# think of a 2D example where we can visualize the local modes, e.g., 
# with two highly correlated predictors. See if the multiple-seed algorithm 
# is about to explore all the local modes (inspiration in simulation of Rockova et al. papers).
