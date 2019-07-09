mice.impute.2l.2stage.norm.intern <-
function(y, ry, x,type,method_est="mm",method_draw="param",incluster=TRUE,kk=5)
{

  
  neighbor<-function(kk,incluster,clust,betalist,xx,ry,yy){
    donneur<-rep(NA,sum(!ry))
    names(donneur)<-seq(nrow(xx))[!ry]
    for (iii in seq(nrow(xx))[!ry]){
      iiifit<-tcrossprod(as.matrix(xx[iii,-which(colnames(xx)==clust),drop=FALSE]),betalist[[xx[iii,clust]]])
      otherfit<-do.call(rbind,mapply(FUN=function(zzz,betazzz){tcrossprod(as.matrix(zzz),betazzz)},zzz=split(xx[,-which(colnames(xx)==clust),drop=FALSE],xx[,clust]),betazzz=betalist))
      if(incluster&sum(ry[xx[,clust]==xx[iii,clust]])>=kk){
        otherfit[xx[,clust]!=xx[iii,clust]]<-Inf
      }else if(incluster&sum(ry[xx[,clust]==xx[iii,clust]])<kk){
        warnings(paste(kk,"neighbours are not available in cluster",xx[iii,clust],". Neighbours are drawn from all clusters."))
      }
      otherfit[!ry]<-Inf
      donneur[as.character(iii)]<-sample(order(abs(as.vector(otherfit)-as.vector(iiifit)))[1:kk],size=1)
    }
    return(yy[donneur])
  }
  
  ana1_imput<-function(dta,fixe)
  {
    fit<-try(lm(formula(paste("y~",paste(fixe[-1],collapse="+"))),data=dta),silent = TRUE)
    if(class(fit)[1]=="try-error") {
      fit<-NA
    }else if(sum(is.na(coef(fit))>0)){
      fit<-NA
      warning("some sporadically missing values imputed as systematically missing")
    }else if(fit$df.residual<=4){
      fit<-NA
      warning("some sporadically missing values imputed as systematically missing\n
              because of overfitting")
    }
    return(fit)
  }
  
  ana_coh3<-function(xx)
  {
    return(list(est=lapply(xx,function(y){if(is.na(y)[1]) NA else coef(y)}),
                varest=lapply(xx,function(y){if(is.na(y)[1]) NA else vcov(y)}),
                sig2=lapply(xx,function(y){if(is.na(y)[1]) NA else summary(y)$sigma^2}),
                ddl=lapply(xx,function(y){if(is.na(y)[1]) NA else y$df.residual}),
                noms= names(xx)))
  }
  
  draw_psi<-function(res.mvmeta,method_est){
    if(method_est=="mm"){
      return(res.mvmeta$Psi)
    }
    
    #reml
    Hes<-as.matrix(forceSymmetric(res.mvmeta$hessian))
    qrsh<-qr.solve(-Hes)
    
    
    valprop<-eigen(qrsh)
    if(sum(valprop$values<0)>0){
      warning("The cov matrix is not full rank")
      valprop$values[valprop$values<0]<-0
      temp<-diag(length(valprop$values));diag(temp)<-valprop$values
      qrsh<-valprop$vectors%*%tcrossprod(temp,valprop$vectors)
    }
    
    #draw
    
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
  names(x) <- paste("V",1:ncol(x),sep="")
  type <- c(2, type)
  
  clust<-names(x)[type==(-2)]
  rande<-names(x)[type==2]
  fixe<-names(x)[type>0]
  
  n.class <- length(unique(x[,clust]))
  
  x[,clust] <- factor(x[,clust], labels=1:n.class)
  
  lev<-levels(x[,clust])
  
  X<-as.matrix(x[,fixe])
  Z<-as.matrix(x[,rande])
  xobs<-x[ry,]
  yobs<-y[ry]
  Xobs<-as.matrix(X[ry,])
  Zobs<-as.matrix(Z[ry,])
  
  fit<-list(NULL)
  beta<-matrix(NA,ncol=length(fixe),nrow=length(lev))
  varbeta<-matrix(NA,ncol=length(fixe)^2,nrow=length(lev))
  sigma<-vector("numeric",length(lev))
  
  temp<-data.frame(y,ry,x)
  temp$y[!temp$ry]<-NA
  
  ####################################
  #meta-analyse 
  ###################################
  #Step 1
  titi<-ana_coh3(by(temp,x[,clust],ana1_imput,fixe=fixe))
  titi<-lapply(titi,FUN=function(x){return(x[!is.na(x)])})
  
  #####################################
  #step2 for sigma
  #####################################
  titi2<-log(matrix(unlist(titi$sig2),nrow=length(titi$sig2),byrow=T))
  S<-lapply(titi$ddl,function(x){2/(x-4)})
  restmp2<-try(mvmeta(titi2~1,S=S,method=method_est,control=list(hessian=TRUE)),silent=T)
  
  if(!is.null(attr(restmp2,"class"))){
    if(attr(restmp2,"class")=="try-error"){
      warning("Pb with mvmeta method=reml, run with method=ml")
      
      restmp2<-mvmeta(titi2~1,S=S,method="ml",control=list(hessian=TRUE))
    }
  }
  
  #draw of sigma star
  res_star<-rnorm(1,mean=coef(restmp2),sd=sqrt(vcov(restmp2)))
  psi_star_sigma2<-draw_psi(res.mvmeta = restmp2,method_est=method_est)
  res_star<-c(res_star,psi_star_sigma2)
  if(is.na(res_star[1])){
    res_star<-restmp2$par
  }
  
  #draw of sigma_i star
  
  ss<-res_star[1]
  WWs<-res_star[2]
  WWsinv<-1/WWs
  ssWWsinv<-ss*WWsinv
  sigma2_i_star<-vector(length=length(lev),mode = "list");names(sigma2_i_star)<-lev
  for (jj in lev)
  {
    if((jj %in% unique(xobs[,clust])) & (!is.null(titi$est[[jj]])))
    {
      
      #sporadically missing values
      
      #sigma
      
      ssi<-log(titi$sig2[[jj]])
      VVsi<-2/(titi$ddl[[jj]]-4)
      sss<-1/(WWsinv+1/VVsi)
      mms<-(ssWWsinv+ssi/VVsi)*sss
      sigma2_i_star[[jj]] <-exp(rnorm(1,mms,sd=sqrt(sss)))
    }
    else{
      #systematically missing
      sigma2_i_star[[jj]] <-exp(rnorm(1,res_star[1],sd=sqrt(res_star[2])))
    }
    
  }
  #####################################
  #step2 for beta and psi
  #####################################
  S<-lapply(sigma2_i_star,function(xx,temp,fixe){xx*solve(crossprod(as.matrix(temp[,fixe])))},temp=temp,fixe=fixe)
  S<-S[names(titi$est)]
  
  titi2<-matrix(unlist(titi$est),nrow=length(titi$est),byrow=T,dimnames=list(names(titi$est),NULL))
  restmp<-try(mvmeta(titi2~1,S=S,method=method_est,control=list(hessian=TRUE)),silent=T)
  
  if(!is.null(attr(restmp,"class")))
  {if(attr(restmp,"class")=="try-error")
  {
    warning("Pb with mvmeta method=reml in 2l.norm.me22, run with method=ml")
    restmp<-mvmeta(titi2~1,S=S,control=list(hessian=TRUE),method="ml")
  }
  }
  
  #draw of beta star
  beta_star<-rmvnorm(1,mean=coef(restmp),sigma=vcov(restmp),method="svd")
  
  #draw of psi
  psi_star_beta<-draw_psi(res.mvmeta = restmp,method_est=method_est)
  
  #draw of beta_i and imputation
  
  WW<-ginv(psi_star_beta)
  WWbb<-WW%*%as.vector(beta_star)
  
  if(method_draw=="param"){
  for (jj in lev)
  {
    if((jj %in% unique(xobs[,clust])) & (!is.null(titi$est[[jj]])))
    {
      
      #sporadically missing
      
      #draw from the posterior distribution of the coefficients
      # beta
      
      bbi<-as.vector(titi$est[[jj]])#hat_beta_i
      VVi<-ginv(titi$varest[[jj]])
      ssi<-ginv(WW+VVi)
      mm<-ssi%*%(WWbb+VVi%*%bbi)
      beta_i_star <- rmvnorm(1,mean=mm,sigma=ssi,method="svd")
      
      
      #imputation
      y[!ry & x[,clust]==jj]<-tcrossprod(X[!ry & x[,clust]==jj,],beta_i_star)+rnorm(sum(!ry & x[,clust]==jj)) * sqrt(sigma2_i_star[[jj]])
      
    }
    else{
      #systematically missing
      beta_i_star <- rmvnorm(1,mean=beta_star,sigma=psi_star_beta,method="svd")
      y[!ry & x[,clust]==jj]<-tcrossprod(X[!ry & x[,clust]==jj,],beta_i_star)+
        rnorm(sum(!ry & x[,clust]==jj)) * sqrt(sigma2_i_star[[jj]])
    }
    
  }
  }else if(method_draw=="pmm"){
    betalist<-vector("list", length(lev))
    names(betalist)<-lev
    for (jj in lev)
    {
      if((jj %in% unique(xobs[,clust])) & (!is.null(titi$est[[jj]])))
      {
        
        #sporadically missing
        
        #draw from the posterior distribution of the coefficients
        # beta
        
        bbi<-as.vector(titi$est[[jj]])#hat_beta_i
        VVi<-ginv(titi$varest[[jj]])
        ssi<-ginv(WW+VVi)
        mm<-ssi%*%(WWbb+VVi%*%bbi)
        betalist[[jj]] <- rmvnorm(1,mean=mm,sigma=ssi,method="svd")
        
        
      }
      else{
        #systematically missing
        betalist[[jj]] <- rmvnorm(1,mean=beta_star,sigma=psi_star_beta,method="svd")

      }
      
    }
    y[!ry]<-neighbor(kk=kk,incluster=incluster,clust=clust,betalist=betalist,xx=x,ry=ry,yy=y)
    }
  return(y[!ry])
}
