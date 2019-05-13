library(parallel)

require(echoseq)
require(locus)
require(R2OpenBUGS)

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
    y <- as.matrix(dat_g$phenos)
    X <- as.matrix(dat_g$snps)
    n <- length(y)
    p <- dim(X)[2]
    d <- dim(y)[2]
    a <- as.vector(single_vb_g$list_hyper$a)
    b <- as.vector(single_vb_g$list_hyper$b)
    lambda <- single_vb_g$list_hyper$lambda
    kappa <- as.vector(single_vb_g$list_hyper$kappa)
    nu <- single_vb_g$list_hyper$nu
    eta <- as.vector(single_vb_g$list_hyper$eta)
    eps <- 1e-3
    WINE="/usr/local/bin/wine"
    WINEPATH="/usr/local/bin/winepath"
    OpenBUGS.pgm="/Applications/OpenBUGS323/OpenBUGS.exe"
    
    data <- c("y","X","p","n", "d","eps", "a", "b","lambda","kappa","eta","nu")
    
    bugs.data(data, dir=getwd(), data.file="data.txt")
    inits <- list(list(gamma=single_vb_g$list_init$gam_vb))
    parameters <- c("beta","gamma")
    bugs.sim <- bugs("data.txt", inits=inits, parameters=parameters, model.file="BUGSmodel.txt", n.chains = 1, n.iter=2000,useWINE = TRUE, WINE=WINE, WINEPATH=WINEPATH, OpenBUGS.pgm = OpenBUGS.pgm, debug=TRUE)
  }
}

