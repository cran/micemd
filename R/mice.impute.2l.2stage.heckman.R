mice.impute.2l.2stage.heckman <- function(y, ry, x, wy=NULL, type, pmm=FALSE, ypmm=NULL, meta_method="reml", pred_std=FALSE,...){
  
   if (pred_std == TRUE){
     x <- apply(x, MARGIN = 2, FUN = function(X) (if(length(unique(X))<=2&sum(unique(X))==1){X} else{scale(X)}))
  }

  # 1. Define variables and dataset----
  
  # Rename covariates
  colnames(x) <- paste0("x_", 1:length(colnames(x))) #change the covariates name for avoiding conflicts when y is covariate
  bos_name <- colnames(x)[type ==  2] # names of variables present in both outcome and selection model
  sel_name <- colnames(x)[type == -3] # names of variables in selection model alone
  out_name <- colnames(x)[type == -4] # names of variables in outcome model alone
  
  # # Define y type
  if ((inherits(y, "factor")) & nlevels(y) == 2){
    family <- "binomial"
  }else{
    family <- "gaussian"
  }
  
  # Check if group variable is defined
  if (length(colnames(x)[type == -2]) == 0) {
    message("No group variable has been provided, the Heckman imputation model will be applied globally to the dataset.")
    Grp_est <- 0 # Group indicator 0: Heckman on full dataset, 1: Heckman at cluster level
  } else {
    group_name <- colnames(x)[type == -2]
    names.clust <- as.character(unique(x[, group_name]))
    Grp_est <- 1
  }
  
  #Define position of selection and outcome in coefficient vector
  order <- c(sel_name, bos_name, out_name)
  send <- length(sel_name) + length(bos_name) + 1 # Index where selection equation ends in coefficient vector
  oend <- send + length(bos_name) + length(out_name) + 1 # Index where outcome equation ends in coefficient vector
  
  # Define outcome and selection equation
  out <- as.formula(paste0("y", "~", paste(c(bos_name,out_name), collapse = "+")))
  sel <- as.formula(paste0("ry", "~", paste(c(sel_name, bos_name), collapse = "+")))
  
  # Define data & prediction matrix
  data <- data.frame(ry, y, x[, order])
  X <- data.frame(cbind(Int=rep(1, nrow(x)), x[, order]))
  
  # 2. Step 1: Get the theta estimators for each study ----
  Syst_nest <- Heck_est <- 1
  if (Grp_est == 1) { #Heckman at cluster level
    df_list <- split(data, x[, group_name]) 
    res_total <- suppressWarnings( lapply( df_list, copulaIPD,
                                           sel = sel, out = out, send = send, family = family)) #Calculate parameters for each study
    studytype <- sapply(res_total, function(x) x[[2]]) # Specify missing pattern in studies
    fit_list  <- lapply(res_total, function(x) x[[1]]) # Get the parameter estimates
    coef_list <- lapply(fit_list[studytype != 0], `[[`, c('coefficients')) # Get effect vector for each study
    Vb_list   <- lapply(fit_list[studytype != 0], `[[`, c('Vb')) # Get covariance matrix for each study
    coef_mat_s <- do.call(rbind, coef_list)
    
    varnam <- colnames(coef_mat_s)
    selnam <- varnam[grepl("*_s" , varnam)]
    outnam <- varnam[!varnam %in% c(selnam, "sigma.star", "theta.star")]
    
    # 3. Step 2: Get marginal theta and var(theta) ----
    
    if( length(studytype[studytype != 0]) < 2 ){
      Grp_est <- 0
    }else{
      Heck_mod <- get_marginal(coef_mat_s, Vb_list, selnam, outnam, meta_method)
      Heck_est <- Heck_mod$Mvma_est}
    Sys.nest <-as.numeric(length(studytype[studytype == 0])==0)
    
  } 
  
  if (Grp_est == 0 | Heck_est == 0){ # Heckman on full dataset or Heckman model no estimable
    message("The Heckman model cannot be estimated marginally, so systematically missing groups will be imputed with the Heckman model based on the full dataset.")
    Heck_mod <- copulaIPD( data = data, sel = sel, out = out, family = family, send = send)
    
    if(Heck_mod[[2]] == 0 &(Grp_est==0|Syst_nest==0)){
      stop("There is insufficient information to impute the Heckman model at the marginal or study level.")
    }
  }
  
  
  # 4. Get theta_k and var(theta_k) from full conditional distribution ----
  if (Grp_est == 1) { # Applies imputation at cluster level
    for (i in names.clust) { #Loop across studies
      
      if (studytype[[i]] == 0) { #systematically missing
        star <- star_systematic( Heck_mod, send, oend, family)
        
      } else { # sporadically missing
        star <- star_sporadic( Heck_mod,
                               coef_list_i = coef_list[[i]],
                               Vb_list_i = Vb_list[[i]], selnam, outnam,family)
      }
      
      Xm <- X[!ry & x[, group_name] == as.numeric(i),]
      
      if (nrow(Xm) != 0) { #Cluster with at least one missing value in outcome equation
        y.star <- gen_y_star( Xm= Xm, sel_name = sel_name, bos_name = bos_name,
                              out_name =out_name, beta_s_star = star$beta_s_star,
                              beta_o_star = star$beta_o_star, sigma_star = star$sigma_star,
                              rho_star = star$rho_star, pmm = pmm, ypmm=ypmm, y = y, ry = ry)
        
        y[!ry & x[, group_name] == as.numeric(i)] <- y.star
      }
    }
  }else{ # Applies imputation on full dataset
    
    star<-star_systematic(Heck_mod, send, oend, family)
    Xm <- X[!ry,]
    
    if (nrow(Xm) != 0) { 
      y.star <- gen_y_star( Xm= Xm, sel_name = sel_name, bos_name = bos_name,
                            out_name =out_name, beta_s_star = star$beta_s_star,
                            beta_o_star = star$beta_o_star, sigma_star = star$sigma_star,
                            rho_star = star$rho_star, pmm = pmm, ypmm=ypmm, y = y, ry = ry)
      y[!ry] <- y.star
    }
  }
  
  return(y[!ry])
}