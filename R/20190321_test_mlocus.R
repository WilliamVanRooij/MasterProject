library(parallel)
require(echoseq)
require(locus)

n <- 200; p <- 500; p0 <- 5; d <- 1; d0 <- 1

cor_type <- "autocorrelated"; 

vec_rho <- runif(floor(p/20), min = 0.95, max = 0.99)

nb_cpus <- 4; # Unfortunately, the number of cores available on Windows is only 1. But when computing on another computer, one can change the number of cpus used.

seed <- 123; set.seed(seed);

list_snps <- generate_snps(n, p, cor_type, vec_rho, n_cpus = nb_cpus,
                           user_seed = seed)

list_phenos <- generate_phenos(n, d, user_seed = seed)

ind_p0 <- sample(1:p, p0)
dat_g <- generate_dependence(list_snps, list_phenos, ind_d0 = sample(1:d, d0),
                             ind_p0 = ind_p0, vec_prob_sh = 0.1,
                             family = "gaussian", max_tot_pve = 0.5, 
                             block_phenos = TRUE, user_seed = seed)

# vb_g <- locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0, link = "identity", user_seed = seed)


user_seed <- sample(1:1e3, 100)#c(123,234,345,456,567,678,789,890,901,012) 

params <- list(Y = dat_g$phenos, X = dat_g$snps, p0_av = 100, link = "identity");

mlocus <- function(fseed) {
  vb_g <- locus(Y = params$Y, X=params$X, p0_av = params$p0_av, link = params$link, user_seed = fseed)
  return(vb_g$gam_vb[,1]) # if d>1 and ind_d0 sampled randomly, this response may be inactive (i.e. no SNP associated with it, so there will be nothing to see...)
}

mac = (Sys.info()['sysname'] != "Windows")

if(!mac) {
  cores <- detectCores()

  cl <- makeCluster(cores)

  mean(unlist(clusterApply(cl=cl,x=user_seed,fun=mlocus))) # are you sure you still get a vector here? maybe see below?

}

# MAC
if(mac) {
  cores <- detectCores()
  
  m_vb_g <- mclapply(user_seed, mlocus, mc.cores = cores) # = n_cpus instead of 1? or = cores if using detectCores?
  
  out <- Reduce('+',  m_vb_g) / length(user_seed) # a bit more compact and no "hard coded" numbers
  plot(out, main='Probabilities of link between a phenotype and SNPs',type='h',lwd=2,lend=1, ylim = c(0,1))
  points(ind_p0, out[ind_p0], col = "red")
  # sum <- 0*(1:p)
  
}
id <- 1
plot(m_vb_g[[id]], main='Probabilities of link between a phenotype and SNPs',type='h',lwd=2,lend=1, ylim = c(0,1))
points(ind_p0, out[ind_p0], col = "red")
# Next steps:
# compare with results from a single seed
# assess the performance with ROC curves, e.g., using the ROCR package
# think of a 2D example where we can visualize the local modes, e.g., 
# with two highly correlated predictors. See if the multiple-seed algorithm 
# is about to explore all the local modes (inspiration in simulation of Rockova et al. papers).
