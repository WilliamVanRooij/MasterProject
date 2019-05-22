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

n <- 100; 
p <- 2; p0 <- 1; 
d <- 1; d0 <- 1


iter <- 1


for (k in sample(1:1e3,iter)){
  
  set.seed(k)  
  
  seed <-  k;
  
  cor_type <- "autocorrelated"; 
  
  #vec_rho <- runif(floor(p/10), min = 0.98, max = 0.99)
  vec_rho <- c(0.98)
  
  nb_cpus <- 4;
  
  ind_d0 <-  sample(1:d, d0)
  
  ind_p0 <- c(1)
  
  p0_av <- 1.5
  
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
  m_vb_g <- mclapply(user_seed, mlocus, mc.cores = nb_cpus)
  
  if(TRUE) {
    
   
    out <- 0
    lb_exp <- 0
    
    elbo <- NULL
    mu_m_locus_1 <- NULL
    mu_m_locus_2 <- NULL
    gam_m_locus_1 <- NULL
    gam_m_locus_2 <- NULL
    sig_m_locus <- NULL
    m_mix_mean_1 <- NULL
    m_mix_mean_2 <- NULL
    m_mix_sig_1 <- NULL
    m_mix_sig_2 <- NULL
    
    
    if(TRUE){
      for(i in c(1:length(user_seed))) {
        out <- out + m_vb_g[[i]]$locus$gam_vb*exp(m_vb_g[[i]]$locus$lb_opt)
        lb_exp <- lb_exp + exp(m_vb_g[[i]]$locus$lb_opt)
        elbo <- c(elbo, exp(m_vb_g[[i]]$locus$lb_opt))
        mu_m_locus_1 <- c(mu_m_locus_1,m_vb_g[[i]]$mu[1],0)
        mu_m_locus_2 <- c(mu_m_locus_2,m_vb_g[[i]]$mu[2],0)
        sig_m_locus <- c(sig_m_locus, m_vb_g[[i]]$sig,0.001)
        gam_m_locus_1 <- c(gam_m_locus_1, exp(m_vb_g[[i]]$locus$lb_opt)*m_vb_g[[i]]$locus$gam_vb[1],exp(m_vb_g[[i]]$locus$lb_opt)*(1-m_vb_g[[i]]$locus$gam_vb[1]))
        gam_m_locus_2 <- c(gam_m_locus_2, exp(m_vb_g[[i]]$locus$lb_opt)*m_vb_g[[i]]$locus$gam_vb[2],exp(m_vb_g[[i]]$locus$lb_opt)*(1-m_vb_g[[i]]$locus$gam_vb[2]))

      }
    }
    
    out <- out / lb_exp
    
  }
  
  single_vb_g <-locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0_av, link = "identity", user_seed = seed, verbose = FALSE, save_hyper=TRUE, save_init = TRUE,full_output=TRUE)
  
  
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
    eps <- 1e4
    
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
    
    if(FALSE){
      mix1 <- norMix(c(single_vb_g$mu_beta_vb[1],0),sigma = c(single_vb_g$sig2_beta_vb,0.01),w=c(single_vb_g$gam_vb[1],1-single_vb_g$gam_vb[1]))
      mix2 <- norMix(c(single_vb_g$mu_beta_vb[2],0),sigma = c(single_vb_g$sig2_beta_vb,0.01),w=c(single_vb_g$gam_vb[2],1-single_vb_g$gam_vb[2]))
      
      joint.density.plot(rnorMix(1e5,mix1),rnorMix(1e5,mix2))
      
      joint.density.plot(out_coda[,"beta[1]"][[1]], out_coda[,"beta[2]"][[1]], Title="Joint Density Plot", contour=T, color=F)
      
    }
  
    if(TRUE){
      m_mix1 <- norMix(mu=mu_m_locus_1,sigma=sig_m_locus,w=1/lb_exp*gam_m_locus_1)
      m_mix2 <- norMix(mu=mu_m_locus_2,sigma=sig_m_locus,w=1/lb_exp*gam_m_locus_2)
      
      #plot(m_mix1)
      
      rmix1 <- rnorMix(1e5,m_mix1)
      rmix2 <- rnorMix(1e5,m_mix2)
      }
    if(FALSE){
      joint.density.plot(rmix1,rmix2,contour=T,color=F)
      
      kill <- kde2d(rmix1,rmix2)
      image(kill)
      #contour(kill,levels=c(0.1,0.3,0.5,0.7,0.9,1.1,1.3),method="simple",xlim=c(-2,2),ylim=c(-2,2))
      contour(kill,method="simple",levels=c(1e-5,1e-3,1e-1,1))
      plot_ly(x=rmix1, y=rmix2)
      plot_ly(x=rmix1, y=rmix2,type="histogram2dcontour")
    }
    if(TRUE){
      
      s <- subplot(
        plot_ly(x = rmix1, type = "histogram"),
        plotly_empty(),
        plot_ly(x = rmix1, y = rmix2, type = "histogram2dcontour"),
        plot_ly(y = rmix2, type = "histogram"),
        nrows = 2, heights = c(0.2, 0.8), widths = c(0.8, 0.2), margin = 0,
        shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE
      )
      p <- layout(s, showlegend = FALSE)
      p
      
    q}
}




