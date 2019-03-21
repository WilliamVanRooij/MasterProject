library(parallel)
require(echoseq)
require(locus)


i <- 1;

n <- 500; p <- 5000; p0 <- 200; d <- 500; d0 <- 400

gam_vb <- 0*c(1:p0);

cor_type <- "autocorrelated"; vec_rho <- runif(100, min = 0.25, max = 0.95)

nb_cpus <- 4; # Unfortunately, the number of cores available on Windows is only 1. But when computing on another computer, one can change the number of cpus used.

seed <- 123; set.seed(seed);

list_snps <- generate_snps(n, p, cor_type, vec_rho, n_cpus = nb_cpus,
                           user_seed = seed)

list_phenos <- generate_phenos(n, d, cor_type = cor_type, vec_rho = vec_rho,
                               n_cpus = nb_cpus, user_seed = seed)

dat_g <- generate_dependence(list_snps, list_phenos, ind_d0 = sample(1:d, d0),
                             ind_p0 = sample(1:p, p0), vec_prob_sh = 0.05,
                             family = "gaussian", max_tot_pve = 0.5, block_phenos = TRUE, user_seed = seed)




user_seed <- cbind(123,234,345,456,567,678,789,890,901,012)

params <- list(Y = dat_g$phenos, X = dat_g$snps, p0_av = p0, link = "identity");

mlocus <- function(seed) {
  vb_g <- locus(Y = params$Y, X=params$X, p0_av = params$p0_av, link = params$link, user_seed = seed)
  return(vb_g$gam_vb[,1])
}


## WINDOWS

if(false) {
  cores <- detectCores()

cl <- makeCluster(cores)

mean(unlist(clusterApply(cl=cl,x=user_seed,fun=mlocus)))

}

# MAC
if(true) {
  
  mclapply(params, locus, mc.cores = nb_cpus)
  
  gam_vb = gam_vb + vb_g$gam_vb[,1]
  
  
  
}


plot(gam_vb)

colPal <- colorRampPalette(c('red','black'))

couleurs <- colPal(500)[gam_vb > 0.5]

plot(gam_vb, col = couleurs)