library(parallel)

require(echoseq)
require(locus)
require(R2OpenBUGS)

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
d <- 2; d0 <- 1


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
  
  # ================================================================================================================================================ #
  #                                                                                                                                                  #
  #                                                                  FUNCTIONS                                                                       #
  #                                                                                                                                                  #
  # ================================================================================================================================================ #
  
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
  
  if(FALSE){
    BUGSmodel <- function(){
      for( j in 1:d ){
        
        for( i in 1:n ){
          y[i, j] ~ dnorm(mu[i,j], tau[j])
          mu[i,j] <- inprod( X[i,1:p], beta[1:p,j] )
        }
        tau[j] ~ dgamma( eta[j], kappa[j] )
        
        tau.beta[j] <- sigma2.inv*tau[j]
        for( k in 1:p ){
          
          beta[k, j] ~ dnorm( 0, prec.beta[k, j] )
          prec.beta[k, j] <- (1-gamma[k, j])*1/eps + gamma[k, j]*tau.beta[j]
          
          gamma[k,j] ~ dbern( omega[k])
        }
      }
      for( k in 1:p ){
        omega[k] ~ dbeta(a[k], b[k])# <- 0.5# ~ dbeta( a[k], b[k] )
      }
      sigma2.inv ~ dgamma(lambda, nu)
    }
    
    write.model(BUGSmodel, "BUGSmodel.txt")
    model.file1 = paste(getwd(),"BUGSmodel.txt", sep="/")
    file.show("BUGSmodel.txt")
    
  }
  
  if(TRUE){
    y <- dat_g$phenos
    X <- dat_g$snps
    n <- as.numeric(nrow(y))
    p <- as.numeric(ncol(X))
    d <- as.numeric(ncol(y))
    a <- single_vb_g$list_hyper$a
    b <- single_vb_g$list_hyper$b
    lambda <- single_vb_g$list_hyper$lambda
    kappa <- single_vb_g$list_hyper$kappa
    nu <- single_vb_g$list_hyper$nu
    eta <- single_vb_g$list_hyper$eta
    eps <- 1e-4
    
    WINE="/usr/local/bin/wine"
    WINEPATH="/usr/local/bin/winepath"
    OpenBUGS.pgm="/Applications/OpenBUGS323/OpenBUGS.exe"
    
    data <- list(y=y,X=X,d=d,p=p,n=n, eps=eps, a=a, b=b,lambda=lambda,kappa=kappa,eta=eta,nu=nu)
    
    inits <- list(list(gamma=apply(single_vb_g$list_init$gam_vb, 2, function(cc) as.numeric(cc>0.5)), beta=single_vb_g$list_init$mu_beta_vb*single_vb_g$list_init$gam_vb, tau=single_vb_g$list_init$tau_vb, sigma2.inv=single_vb_g$list_init$sig2_inv_vb, omega=rbeta(p,single_vb_g$list_hyper$a,single_vb_g$list_hyper$b)))
    parameters <- c("beta","gamma","tau","sigma2.inv","omega")
    
    #bugs.sim <- bugs(data = data, inits=inits, parameters.to.save = parameters, model.file="BUGSmodel.txt", n.chains = 1, n.iter=2000,useWINE = TRUE, WINE=WINE, WINEPATH=WINEPATH, OpenBUGS.pgm = OpenBUGS.pgm, debug=TRUE)
    
    bugs.sim <- R2OpenBUGS::bugs(data = data, inits=inits, parameters.to.save = parameters, 
                                 model.file="BUGSmodel.txt", n.chains = 1, n.iter=50, n.burnin = 10,DIC=TRUE, 
                                 codaPkg = TRUE, working.directory = getwd(), save.history = TRUE, bugs.seed=1,
                                 useWINE = TRUE, WINE=WINE, WINEPATH=WINEPATH, OpenBUGS.pgm = OpenBUGS.pgm, debug=TRUE)
  }
  
  if(FALSE){
    if(FALSE){
      data(schools)
      
      nummodel <- function(){
        for (j in 1:J){
          y[j] ~ dnorm (theta[j], tau.y[j])
          theta[j] ~ dnorm (mu.theta, tau.theta)
          tau.y[j] <- pow(sigma.y[j], -2)}
        mu.theta ~ dnorm (0.0, 1.0E-6)
        tau.theta <- pow(sigma.theta, -2)
        sigma.theta ~ dunif (0, 1000)
      }
      write.model(nummodel, "nummodel.txt")
      model.file1 = paste(getwd(),"nummodel.txt", sep="/")
      file.show("nummodel.txt")
      
      J <- nrow(schools)
      y <- schools$estimate
      sigma.y <- schools$sd
      data <- list ("J", "y", "sigma.y")
    }
    inits <- function(){
      list(theta = rnorm(J, 0, 100), mu.theta = rnorm(1, 0, 100), sigma.theta = runif(1, 0, 100))}
    
    WINE="/usr/local/bin/wine"
    WINEPATH="/usr/local/bin/winepath"
    OpenBUGS.pgm="/Applications/OpenBUGS323/OpenBUGS.exe"
    
    parameters = c("theta", "mu.theta", "sigma.theta")
    
    schools.sim <- bugs(data, inits, model.file = model.file1, parameters=parameters,n.chains = 3, n.iter = 1000, OpenBUGS.pgm=OpenBUGS.pgm, WINE=WINE, WINEPATH=WINEPATH,useWINE=T,debug=T)
    
    print(schools.sim)
    
  }
}

