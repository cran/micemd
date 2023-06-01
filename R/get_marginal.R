get_marginal <- function(coef_mat_s, Vb_list, selnam, outnam, meta_method ){
  
  beta_s = beta_o = rho_t = sigma_t = NA
  Mvma_est <- 1
  
  # separate the set of parameters in beta_out, beta_s and rho
  
  beta_o <- draw_theta_psi_mar( coef_mat_s = coef_mat_s,
                                Vb_list = Vb_list,
                                vnames = outnam,
                                meta_method = meta_method,
                                Mvma_est = Mvma_est)
  beta_s <- draw_theta_psi_mar( coef_mat_s = coef_mat_s,
                                Vb_list = Vb_list,
                                vnames = selnam,
                                meta_method = meta_method,
                                Mvma_est = beta_o[[3]])
  rho_t <- draw_theta_psi_mar( coef_mat_s = coef_mat_s,
                               Vb_list = Vb_list,
                               vnames = "theta.star",
                               meta_method = meta_method,
                               Mvma_est = beta_s[[3]]) #copula package calls atanh(rho) as theta.star
  
  if ("sigma.star"%in%colnames(coef_mat_s)) {
    sigma_t <- draw_theta_psi_mar( coef_mat_s = coef_mat_s,
                                   Vb_list = Vb_list,
                                   vnames = "sigma.star",
                                   meta_method = meta_method,
                                   Mvma_est = rho_t[[3]]) #copula package calls log(sigma) as sigma.star
    Mvma_est <- sigma_t[[3]]
  }
  
  return (list ( beta_s = beta_s,
                 beta_o = beta_o,
                 rho_t = rho_t,
                 sigma_t = sigma_t,
                 Mvma_est = Mvma_est))
}
