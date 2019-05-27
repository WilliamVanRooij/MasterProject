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

n <- 1000; 
p <- 2; p0 <- 1; 
d <- 1; d0 <- 1;

p0_av = 1.7


seed = 123

X <- generate_snps(n, p, cor_type = "equicorrelated", vec_rho=0.9, vec_maf=c(0.4, 0.4), user_seed = seed)$snps
tau <- rgamma(d,shape=1, scale=2)
sigma2.inv <- rgamma(d, shape=1, scale=2)
beta_init <- c(0,rnorm(1,mean=0,sd=sqrt(1/sigma2.inv*tau)))
y <- matrix(rnorm(n, mean=X %*% beta_init, sd=1/sqrt(tau)),nrow=n)
  
  
  mlocus <- function(fseed) {
    
    tau_vb_init <- c(1 / var(y))
    
    sig2_beta_vb_init <- 1 / rgamma(d, shape = 2, rate = 1 / tau_vb_init)
    
    gam_vb_init <-  matrix(rbeta(p * d, shape1 = 1, shape2 = d * (p - p0_av)/p0_av), nrow = p)
    
    
    mu_beta_vb_init <- matrix(rnorm(n = p * d , mean=0, sd = 1), nrow = p)
    
    list_init0 <-  set_init(d,p, gam_vb = gam_vb_init, mu_beta_vb = mu_beta_vb_init, 
                            sig2_beta_vb = sig2_beta_vb_init, tau_vb = tau_vb_init)
    
    vb_g <- locus(Y = y, X=X, p0_av = p0_av, link = "identity", list_init = list_init0, full_output = TRUE)
    
    return(list(locus = vb_g, gam_init = gam_vb_init,mu = mu_beta_vb_init, sig = sig2_beta_vb_init))
  }
  
  user_seed <- sample(1:1e3, 100)
  
  # MAC
  m_vb_g <- mclapply(user_seed, mlocus, mc.cores = 1)
  
  if(TRUE) {
    
    
    out <- 0
    lb_exp <- 0
    
    elbo <- NULL
    gam <- NULL
    mu_beta <- NULL
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
        gam <- c(gam, m_vb_g[[i]]$locus$gam_vb)
        mu_beta <- c(mu_beta, m_vb_g[[i]]$locus$mu_beta_vb)
        lb_exp <- lb_exp + exp(m_vb_g[[i]]$locus$lb_opt)
        elbo <- c(elbo, m_vb_g[[i]]$locus$lb_opt)
        mu_m_locus_1 <- c(mu_m_locus_1,m_vb_g[[i]]$locus$mu_beta_vb[1],0)
        mu_m_locus_2 <- c(mu_m_locus_2,m_vb_g[[i]]$locus$mu_beta_vb[2],0)
        sig_m_locus <- c(sig_m_locus, m_vb_g[[i]]$locus$sig2_beta_vb,0.001)
        gam_m_locus_1 <- c(gam_m_locus_1, m_vb_g[[i]]$locus$gam_vb[1],(1-m_vb_g[[i]]$locus$gam_vb[1]))
        gam_m_locus_2 <- c(gam_m_locus_2, m_vb_g[[i]]$locus$gam_vb[2],(1-m_vb_g[[i]]$locus$gam_vb[2]))
        # gam_m_locus_1 <- c(gam_m_locus_1, exp(m_vb_g[[i]]$locus$lb_opt)*m_vb_g[[i]]$locus$gam_vb[1],exp(m_vb_g[[i]]$locus$lb_opt)*(1-m_vb_g[[i]]$locus$gam_vb[1]))
        # gam_m_locus_2 <- c(gam_m_locus_2, exp(m_vb_g[[i]]$locus$lb_opt)*m_vb_g[[i]]$locus$gam_vb[2],exp(m_vb_g[[i]]$locus$lb_opt)*(1-m_vb_g[[i]]$locus$gam_vb[2]))
        
      }
      vec_w_part <- get_p_m_y(elbo)
      vec_w1 <- rep(vec_w_part, each = 2) * gam_m_locus_1
      vec_w2 <- rep(vec_w_part, each = 2) * gam_m_locus_2
    }
    
    out <- out / lb_exp
    
  }
  
  # ======================================================================== #
  #                                                                          #
  #                               SINGLE                                     #
  #                                                                          #
  # ======================================================================== #
  
  s_tau_vb_init <- c(1 / var(y))
  
  s_sig2_beta_vb_init <- 1 / rgamma(d, shape = 2, rate = 1 /s_tau_vb_init)
  
  s_gam_vb_init <-  matrix(rbeta(p * d, shape1 = 1, shape2 = d * (p - p0_av)/p0_av), nrow = p)

  s_mu_beta_vb_init <- matrix(rnorm(n = p * d , mean=0, sd = 1), nrow = p)
  
  s_list_init <-  set_init(d,p, gam_vb = s_gam_vb_init, mu_beta_vb = s_mu_beta_vb_init, 
                          sig2_beta_vb = s_sig2_beta_vb_init, tau_vb = s_tau_vb_init)
  single_vb_g <-locus(Y = as.matrix(y), X=X, p0_av = p0_av, link = "identity", list_init = s_list_init, verbose = FALSE, save_hyper=TRUE, save_init = TRUE,full_output=TRUE)
  
  
  if(TRUE){
    y <- as.vector(y)
    X <- as.matrix(X)
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
    
    # inits <- list(list(gamma=as.numeric(apply(single_vb_g$list_init$gam_vb, 2, function(cc) as.numeric(cc>0.5))), 
    #                    beta=beta_init, tau=as.numeric(single_vb_g$list_init$tau_vb), 
    #                    sigma2.inv=as.numeric(single_vb_g$list_init$sig2_inv_vb), 
    #                    omega=matrix(rbeta(p,single_vb_g$list_hyper$a,single_vb_g$list_hyper$b),nrow = p)))
    
    inits <- list(list(gamma=as.numeric(apply(single_vb_g$list_init$gam_vb, 2, function(cc) as.numeric(cc>0.5))), 
                       beta=single_vb_g$list_init$mu_beta_vb*single_vb_g$list_init$gam_vb, 
                       tau=as.numeric(single_vb_g$list_init$tau_vb), 
                       sigma2.inv= as.numeric(single_vb_g$sig2_inv_vb), 
                       omega=rbeta(p,single_vb_g$list_hyper$a,single_vb_g$list_hyper$b)))
    
    parameters <- c("beta","gamma","tau","sigma2.inv","omega")
    
    out_bugs <- bugs(data = data, inits=inits, parameters.to.save = parameters, 
                     model.file="BUGSmodelBIS.txt", n.chains = 1, n.iter=5000, n.burnin = 2000, 
                     codaPkg = TRUE, working.directory = getwd(), bugs.seed=13,
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
    j <- 1
    mix1 <- norMix(c(m_vb_g[[j]]$locus$mu_beta_vb[1],0),sigma = c(m_vb_g[[j]]$locus$sig2_beta_vb,0.01),w=c(m_vb_g[[j]]$locus$gam_vb[1],1-m_vb_g[[j]]$locus$gam_vb[1]))
    mix2 <- norMix(c(m_vb_g[[j]]$locus$mu_beta_vb[2],0),sigma = c(m_vb_g[[j]]$locus$sig2_beta_vb,0.01),w=c(m_vb_g[[j]]$locus$gam_vb[2],1-m_vb_g[[j]]$locus$gam_vb[2]))
    solo_rmix1 <- rnorMix(1e5,mix1)
    solo_rmix2 <- rnorMix(1e5,mix2)
    
    #m_mix1 <- norMix(mu=mu_m_locus_1,sigma=sig_m_locus,w=rep(1/200, 200))#1/lb_exp*gam_m_locus_1)
    m_mix1 <- norMix(mu=mu_m_locus_1,sigma=sig_m_locus,w=vec_w1)
    #m_mix2 <- norMix(mu=mu_m_locus_2,sigma=sig_m_locus,w=rep(1/200, 200))#=1/lb_exp*gam_m_locus_2)
    m_mix2 <- norMix(mu=mu_m_locus_2,sigma=sig_m_locus,w=vec_w2)
    
    rmix1 <- rnorMix(1e5,m_mix1)
    rmix2 <- rnorMix(1e5,m_mix2)
  
    d1_solo <- density(solo_rmix1)
    d1_m <- density(rmix1)
    
    d2_solo <- density(solo_rmix2)
    d2_m <- density(rmix2)
    
    plot(d1_solo,xlim=c(-1,0.5),col="blue", main="beta_1 = 0", ylim = c(0, 50)) # ylim = c(0, max(d1_solo, d1_m))) 
    lines(d1_m, lty=2)
    abline(v = beta_init[1],col="red",lty=1)
    
    plot(d2_solo,col="blue", main="beta_2 =/= 0", ylim = c(0, 100)) # ylim = c(0, max(d2_solo, d2_m))) 
    lines(d2_m, lty=2)
    abline(v = beta_init[2],col="red",lty=1)
    
    #joint.density.plot(solo_rmix1,solo_rmix2)
    #joint.density.plot(out_coda[,"beta[1]"][[1]], out_coda[,"beta[2]"][[1]], Title="Joint Density Plot", contour=T, color=F)
    
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
  if(FALSE){
    
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
    
    t <- subplot(
      plot_ly(x = out_coda[,"beta[1]"][[1]], type = "histogram"),
      plotly_empty(),
      plot_ly(x = out_coda[,"beta[1]"][[1]], y = out_coda[,"beta[2]"][[1]], type = "histogram2dcontour"),
      plot_ly(y = out_coda[,"beta[2]"][[1]], type = "histogram"),
      nrows = 2, heights = c(0.2, 0.8), widths = c(0.8, 0.2), margin = 0,
      shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE
    )  
    q <- layout(t,showlegend=F)
    q
  }

if(TRUE){
  par(mfrow=c(1,2))
  hist(out_coda[,"beta[1]"][[1]],breaks=50,freq=F,ylim=c(0,50),col="lightgrey",main=expression(paste("Estimations for ", beta[1])))
  lines(d1_solo,lwd=2,lty=2,col="blue")
  lines(d1_m,lwd=2,lty=2,col="green")
  abline(v=beta_init[1],col="red",lwd=2)
  
  hist(out_coda[,"beta[2]"][[1]],breaks=50,freq=F,col="lightgrey",xlim=c(0.8,1.6),ylim=c(1,15),main=expression(paste("Estimations for ", beta[2])))
  lines(d2_solo,lwd=2,lty=2,col="blue")
  lines(d2_m,lwd=2,lty=2,col="green")
  abline(v=beta_init[2],col="red",lwd=2)
  par(mfrow=c(1,1))
  }



