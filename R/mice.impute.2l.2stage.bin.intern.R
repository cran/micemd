mice.impute.2l.2stage.bin.intern <-
function(y, ry, x,type,method_est="mm")
{
  threshold<-1000
  if((class(y)!="factor")){
    stop("y is not binary, you can use 2l.2stage.norm.reml or convert y as a binary variable")
  }else if( (nlevels(y)!=2)){
    stop("y is not binary, you can use 2l.2stage.norm.reml or convert y as a binary variable")
  }
  
  ana1_imput_bin<-function(dta,fixe)
  {
    if(sum(is.na(dta$y))==length(dta$y)){
      fit<-NA
    }else{
      fit<-try(glm(formula(paste("y~",paste(fixe[-1],collapse="+"))),data=dta,family="binomial"),silent = TRUE)
      if(class(fit)[1]=="try-error") {
        fit<-NA
      }else if((sum(is.na(coef(fit))>0))|any(summary(fit)$coef[,"Std. Error"]>threshold)){
        fit<-NA
        warning("some sporadically missing values imputed as systematically missing")
      }
    }
    return(fit)
  }
  
  ana_coh3_bin<-function(xx)
  {
    return(list(est=lapply(xx,function(y){if(is.na(y)[1]) NA else coef(y)}),
                varest=lapply(xx,function(y){if(is.na(y)[1]) NA else vcov(y)}),
                noms= names(xx)))
  }
  
  draw_psi<-function(res.mvmeta,method_est){
    if(method_est=="mm"){
      return(res.mvmeta$Psi)
    }
    
    Hes<-as.matrix(forceSymmetric(res.mvmeta$hessian))
    qrsh<-qr.solve(-Hes)
    valprop<-eigen(qrsh)
    if(sum(valprop$values<0)>0)
    {
      valprop$values[valprop$values<0]<-0
      temp<-diag(length(valprop$values));diag(temp)<-valprop$values
      qrsh<-valprop$vectors%*%tcrossprod(temp,valprop$vectors)
    }
    
    
    rpar<-rmvnorm(1,mean=res.mvmeta$par,sigma=qrsh,method="svd")
    psi <- matrix(0,ncol( res.mvmeta$Psi),ncol(res.mvmeta$Psi))
    psi[lower.tri(psi,diag=TRUE)] <- rpar
    psi <- tcrossprod(psi)
    
    #check
    temp<-svd(psi)$d[1]
    if( (temp>5) & (temp>(10*svd(res.mvmeta$Psi)$d[1]))){
      warning("Psi cannot be drawn, REML estimation is returned")
      return(res.mvmeta$Psi)
    }
    return(psi)
  }
  
  
  x <- cbind.data.frame(rep(1,nrow(x)), x)
  names(x) <- paste("V",1:(ncol(x)),sep="")
  type <- c(2, type)
  
  clust<-names(x)[type==(-2)]
  rande<-names(x)[type==2]
  fixe<-names(x)[type>0]
  
  n.class <- length(unique(x[,clust]))
  x[,clust] <- factor(x[,clust], labels=1:n.class)
  lev<-levels(x[,clust])
  X<-x[,fixe,drop=F]
  Z<-x[,rande,drop=F]
  xobs<-x[ry,]
  yobs<-y[ry]
  Xobs<-X[ry,,drop=F]
  Zobs<-Z[ry,,drop=F]
  
  fit<-list(NULL)
  beta<-matrix(NA,ncol=length(fixe),nrow=length(lev))
  varbeta<-matrix(NA,ncol=length(fixe)^2,nrow=length(lev))
  
  
  temp<-cbind.data.frame(y,ry,x)
  temp$y[!temp$ry]<-NA
  
  
  ####################################
  #meta-analyse 
  ###################################
  #Step 1
  titi<-ana_coh3_bin(by(temp,x[,clust],ana1_imput_bin,fixe=fixe))
  titi<-lapply(titi,FUN=function(x){return(x[!is.na(x)])})
  
  #step2 for beta and psi
  titi2<-matrix(unlist(titi$est),nrow=length(titi$est),byrow=T,dimnames=list(names(titi$est),NULL))
  
  restmp<-try(mvmeta(titi2~1,S=titi$varest,method=method_est,control=list(hessian=TRUE)))
  if(!is.null(attr(restmp,"class")))
  {if(attr(restmp,"class")=="try-error")
  {
    warning("Pb with mvmeta method=reml in 2l.me22.bin, run with method=ml")
    restmp<-try(mvmeta(titi2~1,S=titi$varest,control=list(hessian=TRUE),method="ml"))
    if(class(restmp)=="try-error"){
      warning("2nd stage cannot be performed, missing values are not imputed")
      return(y[!ry])
    }
  }
  }
  #draw of beta
  
  betaest<-try(rmvnorm(1,mean=coef(restmp),sigma=vcov(restmp),method="svd"))
  if(class(betaest)=="try-error"){betaest<-coef(restmp)}
  
  #draw of psi
  varest<-draw_psi(res.mvmeta = restmp,method_est=method_est)
  
  bb<-as.vector(betaest)
  WW<-ginv(varest)
  WWbb<-WW%*%bb
  
  for (jj in lev)
  {
    # if(length(y[!ry & x[,clust]==jj])==0){next()}#no missing values for the current cluster
    if((jj %in% unique(xobs[,clust]))& (!is.null(titi$est[[jj]])))
    {
      
      #sporadically missing
      
      #draw of parameters from their posterior
      
      bbi<-as.vector(titi$est[[jj]])
      VVi<-ginv(titi$varest[[jj]])
      # print(det(titi$varest[[jj]]))
      ss <- as.matrix(forceSymmetric((ginv(WW + VVi))))
      mm<-ss%*%(WWbb+VVi%*%bbi)
      beta.star <- rmvnorm(1,mean=mm,sigma=ss,method="svd")
      
      #imputation
      
      temp<-tcrossprod(as.matrix(X[!ry & x[,clust]==jj,]),beta.star);
      temp<-exp(temp)/(1+exp(temp))
      imp<-rep(levels(y)[1],length(temp));
      imp[runif(length(y[!ry & x[,clust]==jj]))<temp]<-levels(y)[2]
      y[!ry & x[,clust]==jj]<-imp
      
    }
    else
    {
      #systematically missing
      beta.star <- rmvnorm(1,mean=betaest,sigma=varest,method="svd")
      temp<-tcrossprod(as.matrix(X[!ry & x[,clust]==jj,]),beta.star)
      temp<-exp(temp)/(1+exp(temp))
      imp<-rep(levels(y)[1],length(temp));
      imp[runif(length(y[!ry & x[,clust]==jj]))<temp]<-levels(y)[2]
      y[!ry & x[,clust]==jj]<-imp
    }
    
  }
  
  return(y[!ry])
}
