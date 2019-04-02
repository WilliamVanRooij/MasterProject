library(parallel)

require(echoseq)
require(locus)
require(ROCR)

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
p <- 500; p0 <- 15; 
d <- 1; d0 <- 1

# plot(x=NULL,y=NULL,xlim=c(0,1),ylim=c(0,1))

auc = NULL
c_pred <- NULL
c_lab <- NULL
iter <- 10


for (k in c(1:iter)){

set.seed(k)  

seed <-  k;
  
cor_type <- "autocorrelated"; 

vec_rho <- runif(floor(p/20), min = 0.95, max = 0.99)

nb_cpus <- 4;

ind_d0 <-  sample(1:d, d0)

ind_p0 <- sample(1:p, p0)

user_seed <- sample(1:1e3, 100)

vec_prob_sh <-  0.05 # proba that each SNP will be associated with another active phenotype

max_tot_pve <-  0.5 # max proportion of phenotypic variance explained by the active SNPs

list_snps <- generate_snps(n, p, cor_type, vec_rho, n_cpus = nb_cpus,
                           user_seed = seed)

list_phenos <- generate_phenos(n, d, user_seed = seed)

dat_g <- generate_dependence(list_snps, list_phenos, ind_d0 = ind_d0,
                             ind_p0 = ind_p0, vec_prob_sh = vec_prob_sh,
                             family = "gaussian", max_tot_pve = max_tot_pve, 
                             block_phenos = TRUE, user_seed = seed)

# ================================================================================================================================================ #
#                                                                                                                                                  #
#                                                                  FUNCTIONS                                                                       #
#                                                                                                                                                  #
# ================================================================================================================================================ #


mlocus <- function(fseed) {

  tau_vb_init <- 1 / apply(dat_g$phenos, 2, var)
  
  sig2_beta_vb_init <- 1 / rgamma(d, shape = 2, rate = 1 / tau_vb_init)
  
  gam_vb_init <-  matrix(rbeta(p * d, shape1 = 1, shape2 = 4*d-1), nrow = p)
  
  mu_beta_vb_init <- matrix(rnorm(n = p * d,sd = 100), nrow = p)
  
  list_init <-  set_init(d,p, gam_vb = gam_vb_init, mu_beta_vb = mu_beta_vb_init, 
                         sig2_beta_vb = sig2_beta_vb_init, tau_vb = tau_vb_init)
  
  vb_g <- locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = 100, link = "identity", user_seed = fseed, list_init = list_init)
  return(vb_g)
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
  
  for(i in c(1:length(user_seed))) {
    out <- out + m_vb_g[[i]]$gam_vb*exp(m_vb_g[[i]]$lb_opt)
    lb_exp <- lb_exp + exp(m_vb_g[[i]]$lb_opt)
    
  }
  
  out <- out / lb_exp
  
  
  #out <- Reduce('+',  m_vb_g) / length(user_seed) # a bit more compact and no "hard coded" numbers
  if(FALSE) {
  jpeg(paste("multipleProba",k,".jpg",sep=""), width=800, height=600)
  plot(out, main='Probabilities of link between a phenotype and SNPs',type='h',lwd=1,lend=1, ylim = c(0,1))
  single_vb_g <-locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = 100, link = "identity", user_seed = seed, verbose = FALSE)
  points(ind_p0,single_vb_g$gam_vb[ind_p0], col='black', pch=4)
  points(ind_p0, out[ind_p0], col = "red")
  dev.off()
  }
  c_pred <- cbind(c_pred, out)
  c_lab <- cbind(c_lab, c(1:500) %in% ind_p0)
}

pred <- prediction(c_pred, c_lab)

# perf1 <- performance(pred, "auc")
# auc <- append(auc,perf1@y.values)



}


# plot((1:iter), auc, ylim=c(0,1), pch = 19, main = paste("AUC of ",iter," iterations of the algorithm"), xlab = "Iterations", ylab="AUC")


# single_vb_g <-locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = 100, link = "identity", user_seed = seed, verbose = FALSE)
# points(ind_p0,single_vb_g$gam_vb[ind_p0], col='blue', pch=19)

perf1 <- performance(pred, "tpr","fpr")
plot(perf1,avg="vertical",spread.estimate="stderror",spread.scale=2)





# Next steps:
# implement weights in average - (DONE 27.03.19)
# compare with results from a single seed - (DONE 27.03.19)
# assess the performance with ROC curves, e.g., using the ROCR package (DONE - 27.03.19)
# think of a 2D example where we can visualize the local modes, e.g., 
# with two highly correlated predictors. See if the multiple-seed algorithm 
# is about to explore all the local modes (inspiration in simulation of Rockova et al. papers).
