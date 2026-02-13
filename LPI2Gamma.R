##############################################################
# LPI2Gamma PACKAGE CODE 
# Distribution: Gamma (shape=k, scale=theta)
# Formula: Index = (Q_median - L) / (Q3 - Q1)
###############################################################

library(boot)
library(DT)
library(htmltools)

###############################################################
# 1. Internal: Gamma MLE Helper
###############################################################

LPI2Gamma_mle_helper <- function(x){
  
  # Gamma Log-Likelihood
  nll <- function(par){
    k <- par[1]; lam <- par[2] # Shape, Scale
    if(k <= 0 || lam <= 0) return(1e10)
    -sum(dgamma(x, shape=k, scale=lam, log=TRUE))
  }
  
 
  mean_x <- mean(x)
  var_x  <- var(x)
  start_scale <- var_x / mean_x
  start_shape <- mean_x^2 / var_x
  
  fit <- nlm(nll, p=c(start_shape, start_scale), hessian=TRUE)
  k_hat <- fit$estimate[1]
  lam_hat <- fit$estimate[2]
  vcov_mat <- solve(fit$hessian)
  
  list(shape=k_hat, scale=lam_hat, vcov=vcov_mat)
}

###############################################################
# 2. Internal: Point Estimation (MLE)
###############################################################

LPI2Gamma_point_est <- function(x, p1=0.25, p3=0.75, L){
  
  out <- LPI2Gamma_mle_helper(x)
  k <- out$shape; lam <- out$scale
  
  # Gamma Quantiles (qgamma)
  Q1   <- qgamma(p1, shape=k, scale=lam)
  Qmed <- qgamma(0.5, shape=k, scale=lam)
  Q3   <- qgamma(p3, shape=k, scale=lam)
  
  # FORM??L: (Q_median - L) / (Q3 - Q1)
  Val <- (Qmed - L)/(Q3 - Q1)
  
  list(Value=Val, shape=k, scale=lam, vcov=out$vcov)
}

###############################################################
# 3. Internal: Asymptotic CI (Delta Method for Gamma)
###############################################################

LPI2Gamma_asymptotic <- function(x, p1=0.25, p3=0.75, L, alpha=0.05){
  
  out <- LPI2Gamma_point_est(x, p1, p3, L)
  k <- out$shape; lam <- out$scale; V <- out$vcov
  IndexVal <- out$Value
  
  
  Q1   <- qgamma(p1, shape=k, scale=lam)
  Qmed <- qgamma(0.5, shape=k, scale=lam)
  Q3   <- qgamma(p3, shape=k, scale=lam)
  
  Range <- Q3 - Q1
  
 
  dQ_dlambda <- function(Q, lam) Q/lam
  
  dQ_dk_num <- function(p, k, lam){
    eps <- 1e-5
    q_plus  <- qgamma(p, shape=k+eps, scale=lam)
    q_minus <- qgamma(p, shape=k-eps, scale=lam)
    (q_plus - q_minus) / (2*eps)
  }
  
  dQ1_dlam   <- dQ_dlambda(Q1, lam)
  dQmed_dlam <- dQ_dlambda(Qmed, lam)
  dQ3_dlam   <- dQ_dlambda(Q3, lam)
  
  dQ1_dk     <- dQ_dk_num(p1, k, lam)
  dQmed_dk   <- dQ_dk_num(0.5, k, lam)
  dQ3_dk     <- dQ_dk_num(p3, k, lam)
  
  dI_dQmed <- 1 / Range
  dI_dQ1   <- (Qmed - L) / Range^2  
  dI_dQ3   <- -(Qmed - L) / Range^2 
  
  dI_dk   <- dI_dQ1*dQ1_dk + dI_dQmed*dQmed_dk + dI_dQ3*dQ3_dk
  dI_dlam <- dI_dQ1*dQ1_dlam + dI_dQmed*dQmed_dlam + dI_dQ3*dQ3_dlam
  
  grad <- c(dI_dk, dI_dlam)
  
  varI <- t(grad) %*% V %*% grad
  seI <- sqrt(varI)
  
  CI <- c(IndexVal - qnorm(1-alpha/2)*seI,
          IndexVal + qnorm(1-alpha/2)*seI)
  
  list(Value=IndexVal, SE=seI, LCL=CI[1], UCL=CI[2])
}

###############################################################
# 4. Internal: Bootstrap CI
###############################################################

LPI2Gamma_boot_worker <- function(data, indices, p1, p3, L){
  x <- data[indices]
  tryCatch({
    LPI2Gamma_point_est(x, p1, p3, L)$Value
  }, error = function(e) NA)
}

LPI2Gamma_bootstrap <- function(x, R=2000, p1=0.25, p3=0.75, L, alpha=0.05){
  
  boot_stat <- function(data, indices){
    LPI2Gamma_boot_worker(data, indices, p1, p3, L)
  }
  
  b <- boot(data=x, statistic=boot_stat, R=R)
  
 
  ci <- boot.ci(b, conf = 1-alpha, type=c("norm","basic","perc","bca"))
  
  list(
    Value = b$t0,
    norm = ci$normal[2:3],
    basic = ci$basic[4:5],
    perc  = ci$percent[4:5],
    bca   = ci$bca[4:5]
  )
}

###############################################################
# 5. Internal: Nonparametric Method
###############################################################

LPI2Gamma_nonparam <- function(x, p1=0.25, p3=0.75, L, alpha=0.05){
  
  n <- length(x)
  p_med <- 0.5
  
  q1   <- quantile(x, p1, type=8)
  qmed <- quantile(x, p_med, type=8)
  q3   <- quantile(x, p3, type=8)
  
  Range <- q3 - q1
  IndexVal <- (qmed - L) / Range
  
  f1   <- density(x, from=q1, to=q1)$y[1]
  fmed <- density(x, from=qmed, to=qmed)$y[1]
  f3   <- density(x, from=q3, to=q3)$y[1]
  
  v11 <- p1*(1-p1)/(f1^2)
  v22 <- p_med*(1-p_med)/(fmed^2)
  v33 <- p3*(1-p3)/(f3^2)
  
  v12 <- (p1*(1-p_med))/(f1*fmed) 
  v13 <- (p1*(1-p3))/(f1*f3)      
  v23 <- (p_med*(1-p3))/(fmed*f3) 
  
  Sigma <- matrix(c(v11, v12, v13,
                    v12, v22, v23,
                    v13, v23, v33), 3, 3)
  
  dI_dqmed <- 1 / Range
  dI_dq1   <- (qmed - L) / Range^2
  dI_dq3   <- -(qmed - L) / Range^2
  
  grad <- c(dI_dq1, dI_dqmed, dI_dq3)
  
  varI <- t(grad) %*% Sigma %*% grad / n
  seI  <- sqrt(varI)
  
  CI <- c(IndexVal - qnorm(1-alpha/2)*seI,
          IndexVal + qnorm(1-alpha/2)*seI)
  
  list(Value=IndexVal, SE=seI, LCL=CI[1], UCL=CI[2])
}

###############################################################
# 6. OUTPUT HELPER FUNCTIONS (Print & Viewer)
###############################################################

LPI2Gamma_print <- function(summary_table, x=NULL, alpha=0.05){
  
  conf_pct <- (1-alpha)*100
  
  # Extract values
  mleV   <- summary_table$Value[summary_table$Item=="MLE_Value"]
  asymL  <- summary_table$Value[summary_table$Item=="Asymp_LCL"]
  asymU  <- summary_table$Value[summary_table$Item=="Asymp_UCL"]
  
  normL  <- summary_table$Value[summary_table$Item=="Boot_norm_LCL"]
  normU  <- summary_table$Value[summary_table$Item=="Boot_norm_UCL"]
  
  basicL <- summary_table$Value[summary_table$Item=="Boot_basic_LCL"]
  basicU <- summary_table$Value[summary_table$Item=="Boot_basic_UCL"]
  
  percL  <- summary_table$Value[summary_table$Item=="Boot_perc_LCL"]
  percU  <- summary_table$Value[summary_table$Item=="Boot_perc_UCL"]
  
  bcaL   <- summary_table$Value[summary_table$Item=="Boot_bca_LCL"]
  bcaU   <- summary_table$Value[summary_table$Item=="Boot_bca_UCL"]
  
  npV    <- summary_table$Value[summary_table$Item=="Nonpar_Value"]
  npL    <- summary_table$Value[summary_table$Item=="Nonpar_LCL"]
  npU    <- summary_table$Value[summary_table$Item=="Nonpar_UCL"]
  
  # Header
  cat("=====================================\n")
  cat("          LPI2Gamma RESULTS\n")
  cat(sprintf("   Confidence Level: %.0f%%\n", conf_pct))
  cat("     Dist: Gamma (Shape, Scale)\n")
  cat("=====================================\n\n")
  
  # MLE parameters if x given
  if(!is.null(x)){
    mle_par <- LPI2Gamma_mle_helper(x)
    cat(sprintf("   Shape (k_hat)  : %.4f\n", mle_par$shape))
    cat(sprintf("   Scale (lam_hat): %.4f\n\n", mle_par$scale))
  }
  
  # MLE
  cat(sprintf("MLE Value                : %.4f\n\n", mleV))
  
  # Asymptotic CI
  cat(sprintf("Asymptotic CI            : ( %.4f , %.4f )\n\n", asymL, asymU))
  
  # Bootstrap CIs (ALL)
  cat("Bootstrap CIs:\n")
  cat(sprintf("  Normal CI              : ( %.4f , %.4f )\n", normL, normU))
  cat(sprintf("  Basic CI               : ( %.4f , %.4f )\n", basicL, basicU))
  cat(sprintf("  Percentile CI          : ( %.4f , %.4f )\n", percL, percU))
  cat(sprintf("  BCa CI                 : ( %.4f , %.4f )\n\n", bcaL, bcaU))
  
  # Nonparametric
  cat(sprintf("Nonparametric Value      : %.4f\n", npV))
  cat(sprintf("Nonparametric CI         : ( %.4f , %.4f )\n", npL, npU))
  
  cat("=====================================\n")
}

LPI2Gamma_viewer <- function(summary_table, x, alpha, explanation="LPI2Gamma Analysis") {
  
  mle_par <- LPI2Gamma_mle_helper(x)
  conf_level <- (1 - alpha) * 100
  
  df <- data.frame(
    Item = c(
      "Shape (k_hat)", "Scale (lam_hat)",
      "MLE Value", "Asymptotic CI", 
      "Bootstrap Normal", "Bootstrap Basic", "Bootstrap Percent", "Bootstrap BCa",
      "Nonparametric Value", "Nonparametric CI"
    ),
    Value = c(
      sprintf("%.4f", mle_par$shape), sprintf("%.4f", mle_par$scale),
      sprintf("%.4f", summary_table$Value[summary_table$Item=="MLE_Value"]),
      sprintf("( %.4f , %.4f )", summary_table$Value[summary_table$Item=="Asymp_LCL"], summary_table$Value[summary_table$Item=="Asymp_UCL"]),
      sprintf("( %.4f , %.4f )", summary_table$Value[summary_table$Item=="Boot_norm_LCL"], summary_table$Value[summary_table$Item=="Boot_norm_UCL"]),
      sprintf("( %.4f , %.4f )", summary_table$Value[summary_table$Item=="Boot_basic_LCL"], summary_table$Value[summary_table$Item=="Boot_basic_UCL"]),
      sprintf("( %.4f , %.4f )", summary_table$Value[summary_table$Item=="Boot_perc_LCL"], summary_table$Value[summary_table$Item=="Boot_perc_UCL"]),
      sprintf("( %.4f , %.4f )", summary_table$Value[summary_table$Item=="Boot_bca_LCL"], summary_table$Value[summary_table$Item=="Boot_bca_UCL"]),
      sprintf("%.4f", summary_table$Value[summary_table$Item=="Nonpar_Value"]),
      sprintf("( %.4f , %.4f )", summary_table$Value[summary_table$Item=="Nonpar_LCL"], summary_table$Value[summary_table$Item=="Nonpar_UCL"])
    ),
    stringsAsFactors = FALSE
  )
  

  info_text <- paste0(explanation, " (Confidence Level: ", conf_level, "%)")
  
  explanation_html <- tags$div(
    tags$h3("LPI2Gamma Summary", style="margin-bottom:5px;"),
    tags$p(style="font-size:16px; color:#0056b3; margin-top:0px; font-weight:bold;", info_text),
    tags$hr()
  )
  
  css_fix <- tags$style(HTML(".dataTables_wrapper { height: auto !important; overflow-y: hidden !important; } table.dataTable { width: 100% !important; }"))
  
  return(browsable(tagList(css_fix, explanation_html, DT::datatable(df, rownames = FALSE, options = list(pageLength = nrow(df), scrollY = FALSE, paging = FALSE, dom = 't')))))
}

###############################################################
# 7. MAIN FUNCTION: LPI2Gamma 
###############################################################

LPI2Gamma <- function(x, L, p1=0.25, p3=0.75, R=1000, alpha=0.05){ 
  
  # Calculate Estimates 
  mle      <- LPI2Gamma_point_est(x, p1, p3, L)
  asym     <- LPI2Gamma_asymptotic(x, p1, p3, L, alpha=alpha)
  bootci   <- LPI2Gamma_bootstrap(x, R=R, p1, p3, L, alpha=alpha)
  nonpar   <- LPI2Gamma_nonparam(x, p1, p3, L, alpha=alpha)
  
  TAB <- data.frame(
    Item = c("MLE_Value",
             "Asymp_LCL","Asymp_UCL",
             "Boot_norm_LCL","Boot_norm_UCL",
             "Boot_basic_LCL","Boot_basic_UCL",
             "Boot_perc_LCL","Boot_perc_UCL",
             "Boot_bca_LCL","Boot_bca_UCL",
             "Nonpar_Value","Nonpar_LCL","Nonpar_UCL"),
    Value = c(
      mle$Value,
      asym$LCL, asym$UCL,
      bootci$norm[1], bootci$norm[2],
      bootci$basic[1], bootci$basic[2],
      bootci$perc[1], bootci$perc[2],
      bootci$bca[1], bootci$bca[2],
      nonpar$Value, nonpar$LCL, nonpar$UCL
    )
  )
  
  # AUTO OUTPUT
  LPI2Gamma_print(TAB, x, alpha)
  print(LPI2Gamma_viewer(TAB, x, alpha))
  
  invisible(TAB)
}

###############################################################
# 8. EXAMPLE RUN
###############################################################

set.seed(NULL)
x <- rgamma(500, shape=2.5, scale=5)
L <- 2 


sonuclar <- LPI2Gamma(x, L, p1=0.25, p3=0.75, R=500, alpha = 0.05)