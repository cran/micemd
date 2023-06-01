gen_y_star <- function(Xm, sel_name, bos_name, out_name, beta_s_star, beta_o_star,
                       sigma_star,rho_star, pmm, ypmm, y, ry) {
  
  XOBO <- data.matrix(Xm[,colnames(Xm) %in% c("Int",bos_name,out_name)]) %*% as.vector(beta_o_star)
  XSBS <- data.matrix(Xm[,colnames(Xm) %in% c("Int",sel_name,bos_name)]) %*% as.vector(beta_s_star)
  
  if (!is.na(sigma_star)) { # normal missing variable
    
    Ratio <- (-stats::dnorm(XSBS) / (stats::pnorm(-XSBS)))
    Ratio[is.na(Ratio) | is.infinite(Ratio)] <- 0.0
    y.star <- XOBO + as.numeric(sigma_star) * as.numeric(rho_star) * Ratio +
      rnorm(nrow(XSBS), 0, sd = sigma_star)
    
    if (pmm == TRUE) {
      if (is.null(ypmm)){
        idx <- mice::matchindex(y[ry == 1], y.star)
        y.star <- y[ry == 1][idx]
      }else{
        idx <- mice::matchindex(ypmm, y.star)
        y.star <- ypmm[idx]   
      }
    }
    
  } else { #binomial missing variable
    
    p.star <- pbivnorm::pbivnorm(as.vector(XOBO),-as.vector(XSBS),
                                 -as.numeric(rho_star)) / stats::pnorm(-XSBS)
    p.star[is.na(p.star) | (is.infinite(p.star) & p.star < 0) |p.star < 0.0 |
             p.star == "0"] <- 0.0
    p.star[p.star > 1.0 | p.star == "1" |(is.infinite(p.star) & p.star > 0)] <- 1.0
    
    y.star <-rep(levels(y)[1],nrow(XOBO))
    y.star[runif(nrow(XOBO)) < p.star]<-levels(y)[2]
    
  }
  
  return(y.star)
}