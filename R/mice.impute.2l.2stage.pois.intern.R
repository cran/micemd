mice.impute.2l.2stage.pois.intern <- function(y, ry, x,type,method_est="mm")
{
  threshold<-1000
  if(inherits(y,"factor")){
    stop("y is a factor, you have to use method 2l.2stage.bin or to convert y as a continuous variable")
  }
  
  ana1_imput_bin<-function(dta,fixe)
  {
    if(sum(is.na(dta$y))==length(dta$y)){
      fit<-NA
    }else{
      fit<-try(glm(formula(paste("y~",paste(fixe[-1],collapse="+"))),data=dta,family="poisson"),silent = TRUE)
      if(inherits(fit,"try-error")) {
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
      temp<-diag(length(valprop$values));diag(temp)<-valprop$values#matrice diag avec valeurs propre sur la diag (gere pb diag(vp) si vp=0 )
      qrsh<-valprop$vectors%*%tcrossprod(temp,valprop$vectors)
    }
    
    
    rpar<-rmvnorm(1,mean=res.mvmeta$par,sigma=qrsh,method="svd")#method=svd plus stable
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
  
  titi<-lapply(titi,FUN=function(xx){return(xx[!is.na(xx)])})

  #step2 for beta and psi
  titi2<-matrix(unlist(titi$est),nrow=length(titi$est),byrow=T,dimnames=list(names(titi$est),NULL))

  restmp<-try(mvmeta(titi2~1,S=titi$varest,method=method_est,control=list(hessian=TRUE)))

  if(!is.null(attr(restmp,"class")))
  {if(attr(restmp,"class")=="try-error")
  {
    warning("Pb with mvmeta method=reml in 2l.me22.bin, run with method=ml")
    restmp<-try(mvmeta(titi2~1,S=titi$varest,control=list(hessian=TRUE),method="ml"))
    if("try-error"%in%class(restmp)){
      warning("2nd stage cannot be performed, missing values are not imputed")
      return(y[!ry])
    }
  }
  }
  #draw of beta
 
  betaest<-try(rmvnorm(1,mean=coef(restmp),sigma=vcov(restmp),method="svd"))
  if("try-error"%in%class(betaest)){betaest<-coef(restmp)}
  
  #draw of psi
  varest<-draw_psi(res.mvmeta = restmp,method_est=method_est)
  
  bb<-as.vector(betaest)
  WW<-ginv(varest)
  WWbb<-WW%*%bb

  for (jj in lev)
  {
    # print(jj)
    if(length(y[!ry & x[,clust]==jj])==0){next()}#no missing values for the current cluster
    
    if((jj %in% unique(xobs[,clust]))& (!is.null(titi$est[[jj]])))
    {
      #sporadically missing
      
      #draw of parameters from their posterior
      
      bbi<-as.vector(titi$est[[jj]])
      VVi<-ginv(titi$varest[[jj]])
      ss<-ginv(WW+VVi)
      mm<-ss%*%(WWbb+VVi%*%bbi)
      beta.star <- rmvnorm(1,mean=mm,sigma=ss,method="svd")
      
      #imputation
      
      temp<-tcrossprod(as.matrix(X[!ry & x[,clust]==jj,]),beta.star);
      temp<-exp(temp)
      
      if(nrow(temp)==1){imp<-rpois(1,temp)}else if(nrow(temp)>1){imp<-sapply(temp,function(xx,n){rpois(n,xx)},n=1)}
      
      #if temp is too large, exp(temp) is extremely large and imputed values cannot be drawn.
      if(sum(is.na(imp))>0){imp[is.na(imp)]<-mean(imp,na.rm=TRUE);warning("Convergence issues have been detected")}
      
      y[!ry & x[,clust]==jj]<-imp
      
    }
    else
    {

      #systematically missing
      beta.star <- rmvnorm(1,mean=betaest,sigma=varest,method="svd")
      temp<-tcrossprod(as.matrix(X[!ry & x[,clust]==jj,]),beta.star)
      temp<-exp(temp)
      imp<-sapply(temp,rpois,n=1)
      if(sum(is.na(imp))>0){imp[is.na(imp)]<-mean(imp,na.rm=TRUE);warning("Convergence issues have been detected")}
      
      y[!ry & x[,clust]==jj]<-imp
    }

  }
  
  return(y[!ry])
}