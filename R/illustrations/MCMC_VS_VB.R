a <- 1

rm(list= ls())

setwd("~/Documents/MasterProject/R/illustrations/")
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


seed <- 321
RNGkind("L'Ecuyer-CMRG")
set.seed(seed)

source("fun_utils.R")

# MCMC
n_iter <- 10000
burn_in <- 5000

n <- 300 
p <- 2
stopifnot(p > 1)
vec_pat <- c(TRUE, sample(c(TRUE, FALSE), p - 1, replace = TRUE)) # c(FALSE, TRUE) # pattern for true beta: FALSE = beta = 0, T: beta != 0
vec_pat
maf <- 0.4
cor_type <- c("autocorrelated", "equicorrelated")[2]
rho <- 0.995 # or 0.98

p0_av <- min(sum(vec_pat), p - 0.5)
n_repl <- 100

bool_anneal <- TRUE
if (bool_anneal) {
  anneal <- c(1, 5, 10)
  s_col <- "mediumseagreen"
  m_col <- "red" 
} else {
  anneal <- NULL
  s_col <- "darkblue" 
  m_col <- "orange" 
}

eps <- 1e-4
n_cpus <- 16


list_data <- generate_data(n, p, vec_pat, cor_type, rho, maf, seed) 

X <- list_data$X
y <- list_data$y
beta_true <- list_data$beta

vec_seed <- sample(1:1e3, n_repl)

m_vb <- mclapply(vec_seed, function(seed) {mlocus(anneal = anneal, fseed = seed)}, mc.cores = n_cpus)

vec_elbo <- sapply(seq_along(vec_seed), function(ss) m_vb[[ss]]$lb_opt)
vec_w_part <- get_p_m_y(vec_elbo) 

mat_weights <- sapply(1:p, function(ii) { sapply(seq_along(vec_seed), function(ss) vec_w_part[ss] * c(m_vb[[ss]]$gam_vb[ii,], 1 - m_vb[[ss]]$gam_vb[ii,])) })
mat_mix_mu_beta <- sapply(1:p, function(ii) { sapply(seq_along(vec_seed), function(ss) c(m_vb[[ss]]$mu_beta_vb[ii, ], 0))})
vec_sig_beta <- as.vector(sapply(seq_along(vec_seed), function(ss) c(sqrt(m_vb[[ss]]$sig2_beta_vb), eps)))


if(TRUE){
  
  # in order to get list_hyper and list_init
  #
  s_list_init <- generate_list_init(p, seed) 
  
  s_vb <-locus(Y = y, X=X, p0_av = p0_av, link = "identity", 
               list_init = s_list_init, verbose = FALSE, 
               save_hyper=TRUE, save_init = TRUE, full_output=TRUE, tol = 1e-6, anneal = anneal)
  
  y_mcmc <- as.vector(y) # to avoid overiding the above y and X
  X_mcmc <- as.matrix(X)
  n <- as.numeric(length(y))
  p <- as.numeric(ncol(X))
  a <- as.vector(s_vb$list_hyper$a)
  b <- as.vector(s_vb$list_hyper$b)
  lambda <- as.numeric(s_vb$list_hyper$lambda)
  kappa <- as.numeric(s_vb$list_hyper$kappa)
  nu <- as.numeric(s_vb$list_hyper$nu)
  eta <- as.numeric(s_vb$list_hyper$eta)
  
  ### UNCOMMENT ON MAC
  WINE="/usr/local/bin/wine"
  WINEPATH="/usr/local/bin/winepath"
  OpenBUGS.pgm="/Applications/OpenBUGS323/OpenBUGS.exe"
  
  
  data <- list(y = y_mcmc, X = X_mcmc, p=p, n=n, eps=eps, a=a, b=b, lambda=lambda, kappa=kappa, eta=eta, nu=nu)
  
  inits <- list(list(gamma=as.numeric(apply(s_vb$list_init$gam_vb, 2, function(cc) as.numeric(cc>0.5))), 
                     beta=s_vb$list_init$mu_beta_vb*s_vb$list_init$gam_vb, 
                     tau=as.numeric(s_vb$list_init$tau_vb), 
                     # sigma2.inv= as.numeric(s_vb$sig2_inv_vb), 
                     sigma2.inv = as.numeric(s_vb$list_init$sig2_inv_vb),
                     omega=rbeta(p,shape1 = s_vb$list_hyper$a,shape2 = s_vb$list_hyper$b)))
  
  parameters <- c("beta","gamma","tau","sigma2.inv","omega")
  
  out_bugs <- bugs(data = data, inits=inits, parameters.to.save = parameters, 
                   model.file="BUGSmodelBIS.txt", n.chains = 1, n.iter=n_iter, n.burnin = burn_in, 
                   codaPkg = TRUE, working.directory = getwd(), bugs.seed = 13,
                   useWINE = TRUE, WINE=WINE, WINEPATH=WINEPATH, OpenBUGS.pgm = OpenBUGS.pgm,   ### UNCOMMENT ON MAC
                   debug = F)
  
  out_coda <- read.bugs(out_bugs)
  
  if(FALSE){
   xyplot(out_coda)
   densityplot(out_coda)
  }
}


ss <- 1

list_s_mix <- lapply(1:p, function(ii) norMix(m = c(m_vb[[ss]]$mu_beta_vb[ii,], 0), sigma = c(sqrt(m_vb[[ss]]$sig2_beta_vb), eps), w = c(m_vb[[ss]]$gam_vb[ii,], 1-m_vb[[ss]]$gam_vb[ii,])))
list_m_mix <- lapply(1:p, function(ii) norMix(mu=mat_mix_mu_beta[, ii],sigma = vec_sig_beta, w = mat_weights[, ii]))

names(list_s_mix) <- names(list_m_mix) <- paste0("Predictor_", 1:p)


if(TRUE){
  
  if (p < 4) {
    par(mfrow=c(1, p))
  } else {
    par(mfrow=c(2, ceiling(p/2)))
  }

  for (ii in 1:p) {
    
    plot_densities(ii, s_col = s_col, m_col = m_col, bool_anneal = bool_anneal) 
    
  }
  
  par(mfrow = c(1, 1))
  cor(X)
}



par(mfrow=c(2, 2))

plot_densities(1, s_col = s_col, m_col = m_col, bool_anneal = bool_anneal, xlim = c(-0.2, 2.5), ylim =  c(0, 8), bool_leg = FALSE) 
plot_densities(2, s_col = s_col, m_col = m_col,bool_anneal = bool_anneal, breaks = 60, xlim = c(-0.6, 0.6), ylim =  c(0, 20)) 
plot_densities(3, s_col = s_col, m_col = m_col, bool_anneal = bool_anneal, breaks = 60, xlim = c(-1.5, 0.5), ylim =  c(0, 7), bool_leg = FALSE)
plot_densities(4, s_col = s_col, m_col = m_col, bool_anneal = bool_anneal, breaks = 60, xlim = c(-0.2, 0.6), ylim =  c(0, 25), bool_leg = FALSE)


if(FALSE){
  
  for (ii in 1:p) {
    
    plot(list_s_mix[[ii]], col="blue", main = paste0("beta[", ii, "]"), xlim = c(-abs(beta_true[ii]) - 0.2, abs(beta_true[ii]) + 0.2), lwd = 3, lty = 3, p.norm = FALSE) #, ylim = c(0, 50)) # ylim = c(0, max(d1_solo, d1_m))) 
    lines(list_m_mix[[ii]], lty = 1, p.norm = FALSE) 
    abline(v = beta_true[ii],col="red",lty=1)
    
  }
  
  # joint.density.plot(list_s_mix[[1]], list_s_mix[[4]])
  # joint.density.plot(out_coda[,"beta[1]"][[1]], out_coda[,"beta[4]"][[1]], 
  #                    Title="Joint Density Plot", contour = T, color = F)
  
}

