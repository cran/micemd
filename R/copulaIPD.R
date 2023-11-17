copulaIPD <- function(data, sel, out, family, send) {
  
  fit_ind <- 0 # None model estimable for the cluster
  # A. Estimate Heckman model
  fit <- try(GJRM::gjrm( formula = list(sel, out),
                         data = data,
                         margins = c("probit", ifelse(family=="binomial","probit","N")),
                         model = "BSS",
                         gamlssfit = TRUE,
                         extra.regI = "sED",
                         parscale = TRUE),
             silent = TRUE)
  
  if (!any(inherits(fit, "try-error"))) {
    # model is estimable
    ev <- eigen(fit$fit$hessian, symmetric = TRUE, only.values = TRUE)$values
    convh <- min(ev) > 0 # convergence based on hessian positive definiteness
    convg <- max(abs(fit$fit$gradient)) < 10 # convergence based on abs max gradient
    
    #MAR indication 
    CIcon <- summary(fit)$ CItheta 
    MNAR_ind <- !(CIcon[[1]]>-0.001&CIcon[[2]]<0.001) # exclusion of cases that blow up variance
    
    if (MNAR_ind){
      fit_ind <- 2
      if (convh & convg) {
        # MNAR estimable
        fit_ind <- 1 # Heckman model estimable for the cluster
      }
    }
  }
  
  if( fit_ind == 2){
    fit <- NULL
    gam1 <- try(mgcv::gam(formula=sel,data = data,family = "binomial", method ="REML"))
    gam2 <- try(mgcv::gam(formula=out,data = data,family = family, method ="REML"))
    fit_ind <- 0
    if(!any(inherits(gam1, "try-error"))&!any(inherits(gam2, "try-error"))){
     coefficients <- c(gam1$coefficients,gam2$coefficients)
      vp <- c(diag(gam1$Vp),diag(gam2$Vp))
      if(all(c(!is.na(coefficients),vp!=0))){
        s     <- ifelse(family != "binomial",1,0) 
        ncol1 <- ncol(gam1$Vp)
        ncol2 <- ncol(gam2$Vp)
        Vb    <- matrix(0,ncol = ncol1+ncol2+1+s, nrow = ncol1+ncol2+1+s)
        Vb[1:ncol1,1:ncol1] <- gam1$Vp
        Vb[(ncol1+1):(ncol1+ncol2),(ncol1+1):(ncol1+ncol2)] <- gam2$Vp
        Vb[nrow(Vb),nrow(Vb)]<-1/nrow(data)
        if (family != "binomial") {
          coefficients <- c(coefficients, sigma.star = log(sqrt(gam2$scale)))
          Vb[(nrow(Vb)-1),(nrow(Vb)-1)] <- gam2$V.sp}
        
        fit$coefficients <- c(coefficients,theta.star=0)
        fit$Vb  <- Vb
        fit_ind <- 2}
    }
    
  }
  
  if (fit_ind != 0) {
    names <- c(paste0(names(fit$coefficients)[1:send], "_s"),
               names(fit$coefficients[(send + 1):length(names(fit$coefficients))]))
    names(fit$coefficients) <- names
    colnames(fit$Vb) <- names
    rownames(fit$Vb) <- names
  }else{
    fit <- NA}  
  
  return(list(fit, fit_ind))
  
}