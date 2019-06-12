library(parallel)

require(echoseq)
require(locus)
require(ROCR)
require(nor1mix)
require(plotly)

# rm(list= ls())

RNGkind("L'Ecuyer-CMRG")
seed <- 166
set.seed(seed);


n <- 300; 
p <- 2; p0 <- 1; 
d <- 1; d0 <- 1


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


mlocus <- function(fseed) {
  
  tau_vb_init <- 1 / apply(dat_g$phenos, 2, var)
  
  sig2_beta_vb_init <- 1 / rgamma(d, shape = 2, rate = 1 / tau_vb_init)
  
  gam_vb_init <-  matrix(rbeta(p * d, shape1 = 1, shape2 = 4*d-1), nrow = p)
  
  mu_beta_vb_init <- matrix(rnorm(n = p * d , mean=0, sd = 1), nrow = p)
  
  list_init0 <-  set_init(d,p, gam_vb = gam_vb_init, mu_beta_vb = mu_beta_vb_init, 
                          sig2_beta_vb = sig2_beta_vb_init, tau_vb = tau_vb_init)
  
  vb_g <- locus(Y = dat_g$phenos, X=dat_g$snps, p0_av = p0_av, link = "identity", user_seed = fseed, list_init = list_init0, save_hyper=TRUE, anneal = NULL, save_init = T)
  return(list(locus = vb_g, beta_init = gam_vb_init*mu_beta_vb_init))
}



iter <- 1

cor_type <- "autocorrelated"; 


for(seed in sample(1:1e3,iter)){
  
  set.seed(seed)
  
  vec_rho <- c(0.99)
  
  nb_cpus <- 4;
  
  ind_d0 <- sample(1:d, d0)
  ind_p0 <- sample(1:p, p0)
  
  p0_av <- 1
  
  vec_maf <- runif(p, 0.46, 0.5)
  #vec_maf <- NULL
  
  vec_prob_sh <-  0.05 # proba that each SNP will be associated with another active phenotype
  
  max_tot_pve <-  0.5 # max proportion of phenotypic variance explained by the active SNPs
  
  list_snps <- generate_snps(n, p, cor_type, vec_rho, n_cpus = nb_cpus,
                             user_seed = seed, vec_maf = vec_maf)
  
  list_phenos <- generate_phenos(n, d,  user_seed = seed)
  
  dat_g <- generate_dependence(list_snps, list_phenos, ind_d0 = ind_d0,
                               ind_p0 = ind_p0, vec_prob_sh = vec_prob_sh,
                               family = "gaussian", max_tot_pve = max_tot_pve,
                               block_phenos = TRUE, user_seed = seed)
  
  
  user_seed <- sample(1:1e3, 100)
  
  if(T){

    m_vb_g <- mclapply(user_seed, mlocus, mc.cores = nb_cpus)
    
    
    elbo <- NULL
    gam <- NULL
    out <- 0
    lb_exp <- 0
    beta_init <- NULL
    beta_final <- NULL
    mu_m_locus_1 <- NULL
    mu_m_locus_2 <- NULL
    gam_m_locus_1 <- NULL
    gam_m_locus_2 <- NULL
    sig_m_locus <- NULL
    
    
      for(i in c(1:length(user_seed))) {
        gam <- rbind(gam, as.vector(m_vb_g[[i]]$locus$gam_vb))
        #lb_exp <- lb_exp + exp(m_vb_g[[i]]$locus$lb_opt)
        elbo <- c(elbo, m_vb_g[[i]]$locus$lb_opt)
        beta_init <- cbind(beta_init, m_vb_g[[i]]$beta_init)
        beta_final <- cbind(beta_final, m_vb_g[[i]]$locus$beta_vb)
        mu_m_locus_1 <- c(mu_m_locus_1,m_vb_g[[i]]$locus$list_init$mu_beta_vb[1],0)
        mu_m_locus_2 <- c(mu_m_locus_2,m_vb_g[[i]]$locus$list_init$mu_beta_vb[2],0)
        sig_m_locus <- c(sig_m_locus, m_vb_g[[i]]$locus$list_init$sig2_beta_vb,0.001)
        gam_m_locus_1 <- c(gam_m_locus_1, m_vb_g[[i]]$locus$gam_vb[1],(1-m_vb_g[[i]]$locus$gam_vb[1]))
        gam_m_locus_2 <- c(gam_m_locus_2, m_vb_g[[i]]$locus$gam_vb[2],(1-m_vb_g[[i]]$locus$gam_vb[2]))
        
      }
      
      vec_w1 <- rep(vec_w_part, each = 2) * gam_m_locus_1
      vec_w2 <- rep(vec_w_part, each = 2) * gam_m_locus_2
      
    
    
    vec_w_part <- get_p_m_y(elbo)
    out <- colSums(sweep(gam, 1, vec_w_part, "*"))
    
    
    m_mix1 <- norMix(mu=mu_m_locus_1,sigma=sig_m_locus,w=vec_w1)
    m_mix2 <- norMix(mu=mu_m_locus_2,sigma=sig_m_locus,w=vec_w2)
    
    rmix1 <- rnorMix(1e5,m_mix1)
    rmix2 <- rnorMix(1e5,m_mix2)
    
    d1_m <- density(rmix1)
    
    d2_m <- density(rmix2)
    
    # plot(d1_m,col="blue", main="beta_1 = 0") # ylim = c(0, max(d1_solo, d1_m))) 
    # abline(v = beta_init[1],col="red",lty=1)
    # 
    # plot(d2_m,col="blue", main="beta_2 =/= 0") # ylim = c(0, max(d2_solo, d2_m))) 
    # abline(v = beta_init[2],col="red",lty=1)
    
    data <- data.frame(x=rmix1, y=rmix2)
    points <- data.frame(x=beta_init[1,], y=beta_init[2,])
    ggplot(data, aes(x=x,y=y)) + geom_bin2d(bins=100) + theme_bw()
    
    ggplot(data, aes(x=x,y=y)) + stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") + geom_point(data=points)
    
    plot(beta_init[1,],beta_init[2,], type = "p")
    plot(beta_final[1,], beta_final[2,], col = "red")
    # 
    #   
    # s <- subplot(
    #   plot_ly(x = rmix1, type = "histogram"),
    #   plotly_empty(),
    #   plot_ly(x = rmix1, y = rmix2, type = "histogram2dcontour"),
    #   plot_ly(y = rmix2, type = "histogram"),
    #   nrows = 2, heights = c(0.2, 0.8), widths = c(0.8, 0.2), margin = 0,
    #   shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE
    # )
    # p <- layout(s, showlegend = FALSE)
    # p
    

  }
}


