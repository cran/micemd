cov_mat_vector <- function(cov_mat, vnames) {
  cov_mat[upper.tri(cov_mat)] <- "Up"
  cov_vec <- as.vector(cov_mat[vnames, vnames])
  cov_vec <- as.numeric(cov_vec[cov_vec != "Up"])
  return(cov_vec)
}