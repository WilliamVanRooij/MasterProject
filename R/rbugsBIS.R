library(parallel)

require(echoseq)
require(locus)
require(R2OpenBUGS)
require(coda)
require(lattice)
require(nor1mix)

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

n <- 100; 
p <- 2; p0 <- 1; 
d <- 1; d0 <- 1


iter <- 1


for (k in sample(1:1e3,iter)){
  
  set.seed(k)  
  
  seed <-  k;
  
  cor_type <- "autocorrelated"; 
  
  #vec_rho <- runif(floor(p/10), min = 0.98, max = 0.99)
  vec_rho <- c(0.9)
  
  nb_cpus <- 4;
  
  ind_d0 <-  sample(1:d, d0)
  
  ind_p0 <- c(1)
  
  p0_av <- 1
  
  user_seed <- sample(1:1e3, 100)
  
  vec_prob_sh <-  0.1 # proba that each SNP will be associated with another active phenotype
  
  max_tot_pve <-  0.5 # max proportion of phenotypic variance explained by the active SNPs
  
  list_snps <- generate_snps(n, p, cor_type, vec_rho, n_cpus = nb_cpus,
                             user_seed = seed)
  
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
    
    vb_g <- locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0_av, link = "identity", user_seed = fseed, list_init = list_init0, full_output = TRUE)
    return(list(locus = vb_g, gam_init = gam_vb_init,mu = mu_beta_vb_init, sig = sig2_beta_vb_init))
  }
  
  
  # MAC
  
  if(TRUE) {
    
    m_vb_g <- mclapply(user_seed, mlocus, mc.cores = nb_cpus)
    
    out <- 0
    lb_exp <- 0
    
    if(TRUE){
      for(i in c(1:length(user_seed))) {
        out <- out + m_vb_g[[i]]$locus$gam_vb*exp(m_vb_g[[i]]$locus$lb_opt)
        lb_exp <- lb_exp + exp(m_vb_g[[i]]$locus$lb_opt)
        
        
      }
    }
    
    out <- out / lb_exp
    
  }
  
  single_vb_g <-locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0_av, link = "identity", user_seed = seed, verbose = FALSE, save_hyper=TRUE, save_init = TRUE)
  
  
  if(TRUE){
    y <- as.vector(dat_g$phenos)
    X <- as.matrix(dat_g$snps)
    n <- as.numeric(length(y))
    p <- as.numeric(ncol(X))
    a <- as.vector(single_vb_g$list_hyper$a)
    b <- as.vector(single_vb_g$list_hyper$b)
    lambda <- as.numeric(single_vb_g$list_hyper$lambda)
    kappa <- as.numeric(single_vb_g$list_hyper$kappa)
    nu <- as.numeric(single_vb_g$list_hyper$nu)
    eta <- as.numeric(single_vb_g$list_hyper$eta)
    eps <- 1e-4
    
    WINE="/usr/local/bin/wine"
    WINEPATH="/usr/local/bin/winepath"
    OpenBUGS.pgm="/Applications/OpenBUGS323/OpenBUGS.exe"
    
    
    data <- list(y=y, X=X, p=p, n=n, eps=eps, a=a, b=b, lambda=lambda, kappa=kappa, eta=eta, nu=nu)
    
    inits <- list(list(gamma=apply(single_vb_g$list_init$gam_vb, 2, function(cc) as.numeric(cc>0.5)), beta=single_vb_g$list_init$mu_beta_vb*single_vb_g$list_init$gam_vb, tau=as.numeric(single_vb_g$list_init$tau_vb), sigma2.inv=as.numeric(single_vb_g$list_init$sig2_inv_vb), omega=matrix(rbeta(p,single_vb_g$list_hyper$a,single_vb_g$list_hyper$b),nrow = p)))
    parameters <- c("beta","gamma","tau","sigma2.inv","omega")
    
    out_bugs <- bugs(data = data, inits=inits, parameters.to.save = parameters, 
                                 model.file="BUGSmodelBIS.txt", n.chains = 1, n.iter=10000, n.burnin = 5000, 
                                 codaPkg = TRUE, working.directory = getwd(), bugs.seed=1,
                                 useWINE = TRUE, WINE=WINE, WINEPATH=WINEPATH, OpenBUGS.pgm = OpenBUGS.pgm, debug=F)
    
    out_coda <- read.bugs(out_bugs)
    xyplot(out_coda)
    densityplot(out_coda)
    
    # acfplot(out_coda)
    
    # gelman.diag(out_coda) # Need at least two chains
    # gelman.plot(out_coda)
    
    out_summary <- summary(out_coda, q=c(0.025,0.975))
    out_summary$statistics
    out_summary$quantiles
    
    }
    
    if(TRUE){
      mix <- norMix(c(single_vb_g$list_init$mu_beta_vb[1],0),sigma = c(single_vb_g$list_init$sig2_beta_vb,0.0001),w=c(single_vb_g$gam_vb[1],1-single_vb_g$gam_vb[1]))
      plot(mix)
    }
}

