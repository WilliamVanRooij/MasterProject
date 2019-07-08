library(parallel)

require(echoseq)
require(locus)
require(ROCR)

# rm(list= ls())

RNGkind("L'Ecuyer-CMRG")
seed <- 166
set.seed(seed);


n <- 300; 
p <- 500; p0 <- 50; 
d <- 1; d0 <- 1


bool_anneal <- T
if(bool_anneal) {
  anneal <- c(1, 2, 10)
} else {
  anneal <- NULL
}

cor_type <- "autocorrelated"; 

sigma_delta <- 1


# ================================================================================================================================================ #
#                                                                                                                                                  #
#                                                                  FUNCTIONS                                                                       #
#                                                                                                                                                  #
# ================================================================================================================================================ #


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
  
  gam_vb_init <-  matrix(rbeta(p * d, shape1 = 1, shape2 = 4*d-1), nrow = p)
  
  mu_beta_vb_init <- matrix(rnorm(n = p * d , mean=0, sd = 1), nrow = p)
  
  list_init0 <-  set_init(d,p, gam_vb = gam_vb_init, mu_beta_vb = mu_beta_vb_init, 
                          sig2_beta_vb = sig2_beta_vb_init, tau_vb = tau_vb_init)
  
  vb_g <- locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0_av, link = "identity", user_seed = fseed, list_init = list_init0, save_hyper=TRUE, anneal = NULL)
  return(vb_g)
}

a_mlocus <- function(fseed) {
  
  tau_vb_init <- 1 / apply(dat_g$phenos, 2, var)
  
  sig2_beta_vb_init <- 1 / rgamma(d, shape = 2, rate = 1 / tau_vb_init)
  
  gam_vb_init <-  matrix(rbeta(p * d, shape1 = 1, shape2 = 4*d-1), nrow = p)
  
  mu_beta_vb_init <- matrix(rnorm(n = p * d , mean=0, sd = 1), nrow = p)
  
  list_init0 <-  set_init(d,p, gam_vb = gam_vb_init, mu_beta_vb = mu_beta_vb_init, 
                          sig2_beta_vb = sig2_beta_vb_init, tau_vb = tau_vb_init)
  
  vb_g <- locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0_av, link = "identity", user_seed = fseed, list_init = list_init0, save_hyper=TRUE, anneal = anneal)
  return(vb_g)
}




iter <- 50
{
  
  c_pred <- NULL
  c_lab <- NULL
  
  c_pred_m <- NULL
  c_lab_m <- NULL

  c_pred_w <- NULL
  c_lab_w<- NULL
  
  single_pred <-  NULL
  single_lab <-  NULL
  
  c_pred_a <- NULL
  c_lab_a <- NULL
  
  single_pred_a <-  NULL
  single_lab_a <-  NULL
}

runtime_s_a <- runtime_s <- runtime_m_a <- runtime_m <- 0


for(seed in sample(1:1e3,iter)){
  
  set.seed(seed)
  
  
  vec_rho <- runif(floor(p/10), min = 0.98, max = 0.99)
  #vec_rho <- c(0.99,0.99)
  
  nb_cpus <- 4;
  
  ind_d0 <-  sample(1:d, d0)
  
  #ind_p0 <- c(3,13,17,23,43)
  #ind_p0 <- sample(1:50, p0)  #Seulement pour les plots de probabilités
  ind_p0 <- sample(1:p, p0)
  
  p0_av <- 50
  
  vec_maf <- runif(p, 0.4, 0.5)
  #vec_maf <- NULL
  
  vec_prob_sh <-  0.05 # proba that each SNP will be associated with another active phenotype
  
  max_tot_pve <-  0.8 # max proportion of phenotypic variance explained by the active SNPs
  
  list_snps <- generate_snps(n, p, cor_type, vec_rho, n_cpus = nb_cpus,
                             user_seed = seed, vec_maf = vec_maf)
  
  list_phenos <- generate_phenos(n, d,  user_seed = seed)
  
  dat_g <- generate_dependence(list_snps, list_phenos, ind_d0 = ind_d0,
                               ind_p0 = ind_p0, vec_prob_sh = vec_prob_sh,
                               family = "gaussian", max_tot_pve = max_tot_pve,
                               block_phenos = TRUE, user_seed = seed)
  
  
  user_seed <- sample(1:1e3, 100)
  
  if(T){
    time0_m_a <- proc.time()
    m_vb_g_a <- mclapply(user_seed, a_mlocus, mc.cores = nb_cpus)
    runtime_m_a <- runtime_m_a + proc.time() - time0_m_a
    time0_m <- proc.time()
    m_vb_g <- mclapply(user_seed, mlocus, mc.cores = nb_cpus)
    runtime_m <- runtime_m + proc.time() - time0_m
    
    elbo <- NULL
    gam <- NULL
    out <- 0
    lb_exp <- 0
    
    elbo_a <- NULL
    gam_a <- NULL
    out_a <- 0
    lb_exp_a <- 0
    
    out_m <- 0
    elbo_w <- NULL
    w <-  NULL
    
    if(TRUE){
      for(i in c(1:length(user_seed))) {
        gam <- rbind(gam, as.vector(m_vb_g[[i]]$gam_vb))
        #lb_exp <- lb_exp + exp(m_vb_g[[i]]$locus$lb_opt)
        elbo <- c(elbo, m_vb_g[[i]]$lb_opt)
      }
    }
    
    if(TRUE){
      for(i in c(1:length(user_seed))) {
        gam_a <- rbind(gam_a, as.vector(m_vb_g_a[[i]]$gam_vb))
        #lb_exp <- lb_exp + exp(m_vb_g[[i]]$locus$lb_opt)
        elbo_a <- c(elbo_a, m_vb_g_a[[i]]$lb_opt)
      }
    }
    
    if(TRUE){    # Equi-weighted averaged
      for(i in c(1:length(user_seed))){
        out_m <- out_m + m_vb_g[[i]]$gam_vb
        w <- c(w,exp(-sigma_delta*(sum(as.numeric(m_vb_g[[i]]$gam_vb > mean(m_vb_g[[i]]$gam_vb)))-p0_av)^2))
        elbo_w <- c(elbo_w, m_vb_g[[i]]$lb_opt*exp(-sigma_delta*(sum(as.numeric(m_vb_g[[i]]$gam_vb > mean(m_vb_g[[i]]$gam_vb)))-p0_av)^2))
      }
      
      out_m <-  out_m/length(user_seed)
    }
    
    vec_w_part <- get_p_m_y(elbo)
    out <- colSums(sweep(gam, 1, vec_w_part, "*"))
    
    vec_w_part_a <- get_p_m_y(elbo_a)
    out_a <- colSums(sweep(gam_a, 1, vec_w_part_a, "*"))
    
    vec_w_part_w <- get_p_m_y(elbo_w)
    out_w <- colSums(sweep(gam, 1, vec_w_part_w, "*"))
    
    
    time0_s_a <- proc.time()
    single_vb_g_a <-locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0_av, link = "identity", user_seed = seed, verbose = FALSE, save_hyper=TRUE, anneal = anneal)
    runtime_s_a <- runtime_s_a + proc.time() - time0_s_a
    time0_s <- proc.time()
    single_vb_g <-locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0_av, link = "identity", user_seed= seed, verbose = FALSE, save_hyper=TRUE, anneal = NULL)
    runtime_s <- runtime_s + proc.time() - time0_s
    
    single_pred <- cbind(single_pred, single_vb_g$gam_vb)
    single_lab <- cbind(single_lab, c(1:500) %in% ind_p0)
    
    c_pred <- cbind(c_pred, out)
    c_lab <- cbind(c_lab, c(1:500) %in% ind_p0)
    
    single_pred_a <- cbind(single_pred_a, single_vb_g_a$gam_vb)
    single_lab_a <- cbind(single_lab_a, c(1:500) %in% ind_p0)
    
    c_pred_a <- cbind(c_pred_a, out_a)
    c_lab_a <- cbind(c_lab_a, c(1:500) %in% ind_p0)
    
    c_pred_m <- cbind(c_pred_m, out_m)
    c_lab_m <- cbind(c_lab_m, c(1:500) %in% ind_p0)
    
    c_pred_w <- cbind(c_pred_w, out_w)
    c_lab_w <- cbind(c_lab_w, c(1:500) %in% ind_p0)
  }
}

if(T){ # ROC CURVES
  pred_m_locus <- prediction(c_pred, c_lab)
  pred_s_locus <- prediction(single_pred, single_lab)
  
  pred_m_locus_m <- prediction(c_pred_m, c_lab_m)
  perf_m_locus_m <- performance(pred_m_locus_m, "tpr","fpr")
  
  pred_m_locus_w <- prediction(c_pred_w, c_lab_w)
  perf_m_locus_w <- performance(pred_m_locus_w, "tpr","fpr")
  
  perf_m_locus <- performance(pred_m_locus, "tpr","fpr")
  perf_s_locus <- performance(pred_s_locus, "tpr","fpr")
  
  pred_m_locus_a <- prediction(c_pred_a, c_lab_a)
  pred_s_locus_a <- prediction(single_pred_a, single_lab_a)
  
  perf_m_locus_a <- performance(pred_m_locus_a, "tpr","fpr")
  perf_s_locus_a <- performance(pred_s_locus_a, "tpr","fpr")
  
  par(pty="s")
  #pdf(paste("ROC_Comp_p0_",p0,"_var_0_",floor(10*max_tot_pve),".pdf",sep=""))
  plot(perf_m_locus,avg="vertical",spread.estimate="stderror",spread.scale=2,col='orange',lwd=2, main=expression(paste("ROC Curves comparison, ",p[0]," = 50, Max Tot. PVE = 0.8")),xlim=c(0,0.2))
  plot(perf_s_locus,avg="vertical",spread.estimate="stderror",spread.scale=2,col='blue', lwd=2, add=T)
  plot(perf_m_locus_a,avg="vertical",spread.estimate="stderror",spread.scale=2,col='red', lwd=2, add=T)
  plot(perf_s_locus_a,avg="vertical",spread.estimate="stderror",spread.scale=2,col='green', lwd=2, add=T)
  plot(perf_m_locus_m, avg="vertical", spred.estimate="stderror",spread.scale=2,col='purple',lwd=2,add=T)
  #plot(perf_m_locus_w, avg="vertical", spred.estimate="stderror",spread.scale=2,col='darkgreen',lwd=2,add=T)
  legend(0.05,1.0, c("LOCUS","Annealed LOCUS", "Averaged LOCUS","Averaged annealed LOCUS","Averaged LOCUS (Equal weights)"), col=c('blue', 'green','orange', 'red','purple'),lwd=2)
  #dev.off()
  
  par(pty="m")
  
}
par(mfrow=c(1,1))
if(F){
  #make_ld_plot(dat_g$snps[,1:50],"r")
  #png("m_locus.png", width=715, height=350 )
  plot(out[1:50],type='h',lwd=10,lend=1,xlab='',col='#a4a4a4', xaxt='n',main ="Probability of association - Multiple LOCUS", ylab='')
  points(ind_p0, out[ind_p0], col = "red",type='h',lend=1,lwd=10)
  #dev.off()
}

if(F){
  #make_ld_plot(dat_g$snps[,1:50],"r")
  #png("s_locus.png", width=715, height=350 )
  plot(single_vb_g$gam_vb[1:50],type='h',lwd=10,lend=1,xlab='',col='#a4a4a4', xaxt='n',main ="Probability of association - Single LOCUS", ylab='')
  points(ind_p0, single_vb_g$gam_vb[ind_p0],col='red', type='h', lend=1,lwd=10)
  #dev.off()
}

if(F){
  #make_ld_plot(dat_g$snps[,1:50],"r")
  #png("m_annealed.png", width=715, height=350 )
  plot(out_a[1:50],type='h',lwd=10,lend=1,xlab='',col='#a4a4a4', xaxt='n',main ="Probability of association - Annealed multiple LOCUS", ylab='')
  points(ind_p0, out_a[ind_p0], col = "red",type='h',lend=1,lwd=10)
  #dev.off()
}

if(F){
  #make_ld_plot(dat_g$snps[,1:50],"r")
  #png("s_annealed.png", width=715, height=350 )
  plot(single_vb_g_a$gam_vb[1:50],type='h',lwd=10,lend=1,xlab='',col='#a4a4a4', xaxt='n',main ="Probability of association - Annealed single LOCUS", ylab='')
  points(ind_p0, single_vb_g_a$gam_vb[ind_p0],col='red', type='h', lend=1,lwd=10)
  #dev.off()
}
if(T){
  #make_ld_plot(dat_g$snps[,1:50],"r")
  #png("m_annealed.png", width=715, height=350 )
  plot(out_m[1:50],type='h',lwd=10,lend=1,xlab='',col='#a4a4a4', xaxt='n',main ="Probability of association - Annealed multiple LOCUS", ylab='')
  points(ind_p0, out_m[ind_p0], col = "red",type='h',lend=1,lwd=10)
  #dev.off()
}
if(F){
  #make_ld_plot(dat_g$snps[,1:50],"r")
  #png("m_annealed.png", width=715, height=350 )
  plot(out_w[1:50],type='h',lwd=10,lend=1,xlab='',col='#a4a4a4', xaxt='n',main ="Probability of association - Meaned multiple LOCUS", ylab='')
  points(ind_p0, out_w[ind_p0], col = "red",type='h',lend=1,lwd=10)
  #dev.off()
}


if(F){
  #pdf("runtimes.pdf")
  plot(c(runtime_s[2]/iter,runtime_s_a[2]/iter,runtime_m[2]/iter,runtime_m_a[2]/iter),type="h",lwd=50,col=c("blue","green","orange","red"),lend=1, main="Running times of the four methods", xaxt="n",xlab="",ylab="Runtimes",ylim=c(0.002,0.028))
  text(x=c(1,2,3,4),y=c(runtime_s[2]/iter+0.0007,runtime_s_a[2]/iter+0.0007,runtime_m[2]/iter+0.0007,runtime_m_a[2]/iter+0.0007),labels=round(c(runtime_s[2]/iter,runtime_s_a[2]/iter,runtime_m[2]/iter,runtime_m_a[2]/iter),3))
  legend(1.1,0.025, c("LOCUS","Annealed LOCUS", "Averaged LOCUS","Averaged annealed LOCUS"), col=c('blue', 'green','orange', 'red'),lwd=4)
  #dev.off()
}
