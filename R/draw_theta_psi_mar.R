draw_theta_psi_mar <- function(coef_mat_s, Vb_list, meta_method, Mvma_est, vnames = NULL) {
  
  theta_star <- NA
  psi_star <- NA
  
  if (is.null(vnames)) { #use all set of parameters
    vnames <- colnames(coef_mat_s)
  }
  
  # Get covariance matrix
  coef_mat_s <- coef_mat_s[, vnames]
  cov_mat_s <- do.call("rbind", lapply(Vb_list, cov_mat_vector, vnames = vnames))
  
  # Apply multivariate random-effects meta-analysis
  mvma <- suppressWarnings(try(mixmeta::mixmeta(coef_mat_s,cov_mat_s, method = meta_method,
                                                control = list(hessian = TRUE)), silent = TRUE))
  
  if (inherits(mvma,"try-error")) { # Use mm instead
    meta_method = "mm"
    mvma <- suppressWarnings(try(mixmeta::mixmeta(coef_mat_s, cov_mat_s, method = meta_method,
                                                  control = list(hessian = TRUE)),silent = TRUE))
    
    if (inherits(mvma,"try-error")) { # MA can not be estimated
      Mvma_est <- 0}
  }
  
  
  if (Mvma_est == 1) {
    # Draw effects theta_star
    theta_star <- MASS::mvrnorm(n = 1, mu = coef(mvma), Sigma = vcov(mvma))
    
    if (meta_method != "mm") {
      # Draw random effect, psi_star
      if (length(vnames) == 1) {
        qrsh <- 1 / mvma$hessian
      } else {
        Hes <- as.matrix(Matrix::forceSymmetric(mvma$hessian))
        qrsh <- as.matrix(Matrix::nearPD(MASS::ginv(-Hes))$mat)
      }
      
      rpar <- mvtnorm::rmvnorm(1, mean = mvma$par, sigma = qrsh, method = "svd")
      
      if (length(vnames) == 1) {
        psi_star <- rpar ^ 2
      } else {
        psi <- matrix(0, ncol(mvma$Psi), ncol(mvma$Psi))
        psi[lower.tri(psi, diag = TRUE)] <- rpar
        psi_star <- Matrix::tcrossprod(psi)
        
      }
      
    } else # meta_method== reml OR ml
      psi_star <- mvma$Psi
  }
  
  colnames(psi_star) <- names(theta_star)
  rownames(psi_star) <- names(theta_star)
  
  return(list(theta_star, psi_star, Mvma_est))
  
}