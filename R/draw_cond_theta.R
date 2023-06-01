draw_cond_theta <- function(theta_mar, theta_k, var_theta_k, vnames) {
  W_m <- MASS::ginv(theta_mar[[2]])
  W_k <- MASS::ginv(var_theta_k[vnames, vnames])
  S <- MASS::ginv(W_m + W_k)
  mu <- S %*% (W_k %*% as.vector(theta_k[vnames]) + W_m %*% as.vector(theta_mar[[1]]))
  theta_star_i <- MASS::mvrnorm(n = 1, mu = as.vector(mu), Sigma = S)
  return(theta_star_i)
}