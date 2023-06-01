star_sporadic <- function(Heck_mod, coef_list_i, Vb_list_i, selnam, outnam, family){
  
  beta_s_star = beta_o_star = sigma_star = rho_star = NA
  cond <- ifelse(is.null(Heck_mod$Mvma_est),0,Heck_mod$Mvma_est)
  
  if (cond == 1){ # Draw study parameters from conditional distribution
    beta_o_star <- draw_cond_theta( theta_mar = Heck_mod$beta_o, theta_k = coef_list_i,
                                    var_theta_k = Vb_list_i, vnames = outnam)
    beta_s_star <- draw_cond_theta( theta_mar = Heck_mod$beta_s,theta_k = coef_list_i,
                                    var_theta_k = Vb_list_i, vnames = selnam)
    
    rho_star <- tanh( draw_cond_theta( theta_mar = Heck_mod$rho_t,
                                       theta_k = coef_list_i,
                                       var_theta_k = Vb_list_i,
                                       vnames = "theta.star")) #copula package calls atanh(rho) as theta.star
    
    if (family == "gaussian") {
      sigma_star <- exp(draw_cond_theta( theta_mar = Heck_mod$sigma_t,
                                         theta_k = coef_list_i,
                                         var_theta_k = Vb_list_i,
                                         vnames = "sigma.star"))#copula package calls log(sigma) as sigma.star
    }
    
    
  } else { # Draw from study parameters from study distribution
    
    beta_o_star <- MASS::mvrnorm(n = 1,mu = coef_list_i[outnam],Sigma = Vb_list_i[outnam, outnam])
    beta_s_star <- MASS::mvrnorm(n = 1,mu = coef_list_i[selnam],Sigma = Vb_list_i[selnam, selnam])
    rho_star <- tanh(MASS::mvrnorm( n = 1, mu = coef_list_i["theta.star"],
                                    Sigma = Vb_list_i["theta.star", "theta.star"]))
    if (family == "gaussian") {
      sigma_star <- exp(MASS::mvrnorm( n = 1, mu = coef_list_i["sigma.star"],
                                       Sigma = Vb_list_i["sigma.star", "sigma.star"]))#copula package calls log(sigma) as sigma.star
    }
    
  }
  
  return (list ( beta_s_star = beta_s_star,
                 beta_o_star = beta_o_star,
                 rho_star = rho_star,
                 sigma_star = sigma_star))
}