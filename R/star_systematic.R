star_systematic <- function(Heck_mod, send, oend, family){
  
  if (is.null(Heck_mod$Mvma_est)){ # From total data model
    star <- mvtnorm::rmvnorm( 1, mean = Heck_mod[[1]]$coefficients,
                              sigma = Heck_mod[[1]]$Vb, method = "svd")
    
    beta_o_star <- star[(send + 1):oend]
    beta_s_star <- star[1:send]
    rho_star <- tanh(star[, "theta.star"])
    sigma_star <- ifelse(family=="gaussian",
                         exp(star[, "sigma.star"]), NA)
    
    
  } else { # From meta model
    
    beta_o_star <- MASS::mvrnorm(n = 1, mu = Heck_mod$beta_o[[1]],Sigma =Heck_mod$beta_o[[2]])
    beta_s_star <- MASS::mvrnorm(n = 1, mu = Heck_mod$beta_s[[1]],Sigma = Heck_mod$beta_s[[2]])
    rho_star <- tanh(MASS::mvrnorm( n = 1, mu = Heck_mod$rho_t[[1]],Sigma = Heck_mod$rho_t[[2]]))
    if (family =="gaussian") {
      sigma_star <- exp(MASS::mvrnorm( n = 1, mu = Heck_mod$sigma_t[[1]],
                                       Sigma = Heck_mod$sigma_t[[2]])) #copula package calls log(sigma) as sigma.star
    } else { #binomial
      sigma_star <- NA}
  }
  
  return (list (beta_s_star = beta_s_star,
                beta_o_star = beta_o_star,
                rho_star = rho_star,
                sigma_star = sigma_star))}