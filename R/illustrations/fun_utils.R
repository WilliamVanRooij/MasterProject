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

generate_data <- function(n, p, vec_pat, cor_type, rho, maf, seed) {
  
  stopifnot(is.logical(vec_pat))
  
  X <- generate_snps(n, p, 
                     cor_type = cor_type, 
                     vec_rho=rho, vec_maf=rep(maf, p), 
                     user_seed = seed)$snps 
  
  tau <- rgamma(1, shape=1, scale=2)
  sigma2.inv <- rgamma(1, shape=2, scale=1)
  beta <- rep(0, p)
  beta[vec_pat] <- rnorm(sum(vec_pat), mean=0, sd=sqrt(1/(sigma2.inv*tau)))
  y <- matrix(rnorm(n, mean=X %*% beta, sd=1/sqrt(tau)), nrow = n, ncol = 1)
  
  list("X" = X, "y" = y, "beta" = beta)
}

generate_list_init <- function(p, seed) {
  
  set.seed(seed)
  
  tau_vb_init <- rgamma(1, shape = 10, rate = 1)
  
  sig2_beta_vb_init <- 1 / rgamma(1, shape = 1, rate = 1)
  
  gam_vb_init <-  matrix(rbeta(p, shape1 = 0.8, shape2 = 1), nrow = p)
  
  mu_beta_vb_init <- matrix(rnorm(n = p, mean=0, sd = 10), nrow = p)
  
  set_init(d = 1, p, gam_vb = gam_vb_init, mu_beta_vb = mu_beta_vb_init, 
           sig2_beta_vb = sig2_beta_vb_init, tau_vb = tau_vb_init)
  
}


mlocus <- function(anneal, fseed) {
  
  list_init <- generate_list_init(p, fseed)
  
  vb <- locus(Y = y, X=X, p0_av = p0_av, link = "identity", #list_hyper = list_hyper,
        list_init = list_init, full_output = TRUE, tol = 1e-6, anneal = anneal, user_seed = fseed)
  
}


plot_densities <- function(ii, s_col, m_col, bool_anneal, breaks = 50, 
                           xlim = c(-abs(beta_true[ii]) - 1.3, abs(beta_true[ii]) + 1.3), ylim =  c(0, 8), bool_leg = TRUE) {
  
  hist(out_coda[, paste0("beta[", ii, "]")][[1]], breaks = breaks, freq = FALSE, col = "gray90", 
       main = paste0("Distribution of beta[", ii, "]"), 
       xlab = paste0("beta[", ii, "]"), xlim = xlim, 
       ylim = ylim)
  lines(list_s_mix[[ii]], lwd = 2.5, lty=1, col = s_col, p.norm=FALSE)
  lines(list_m_mix[[ii]], lwd = 2, lty=1, col = m_col, p.norm = FALSE)
  abline(v=beta_true[ii], col="black", lty = 2, lwd = 1.5)
  
  if (bool_leg) {
    legend("topright", 
           c(paste0(ifelse(bool_anneal, "Annealing ", ""), c("Single LOCUS", "Multiple LOCUS")), "Simulated beta coefficient"), 
           col = c(s_col, m_col, "black"), lty = c(1, 1, 2), lwd = 2, bty = "n")
  }

}

