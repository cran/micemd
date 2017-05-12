mice.impute.2l.glm.norm <-
function(y, ry, x,type,...){
  # the main code
  x<-cbind.data.frame(rep(1,nrow(x)), x)
  names(x) <- paste("V",1:ncol(x),sep="")
  
  type <- c(2, type)
  
  clust <- names(x)[type==(-2)]
  rande <- names(x)[type==2]
  fixe <- names(x)[type>0]
  
  n.class <- length(unique(x[,clust]))
  x[,clust] <- factor(x[,clust], labels=1:n.class)
  lev<-levels(x[,clust])
  X<-x[,fixe,drop=F]
  Z<-x[,rande,drop=F]
  xobs<-x[ry,,drop=F]
  yobs<-y[ry]
  Xobs<-X[ry,,drop=F]
  Zobs<-Z[ry,,drop=F]
  
  if(length(rande)>1){randmodel <- paste("yobs ~ ", paste(fixe[-1], collapse="+"),
                     "+ ( 1 +", paste(rande[-1],collapse="+"), 
                     "|", clust, ")") # [-1] to remove intercept
  }else{
    randmodel <- paste("yobs ~ ", paste(fixe[-1], collapse="+"),
                       "+ ( 1 ", 
                       "|", clust, ")")
    }
  suppressWarnings(fit <- try(lmer(formula(randmodel), 
                                   data = data.frame(yobs,xobs)),silent=F))
  if(!is.null(attr(fit,"class"))){
    if(attr(fit,"class")=="try-error"){
      warning("lmer cannot be run, sorry!")
      return(y[!ry])
    }
  }  
  
  # draw beta*
  beta <- fixef(fit)
  rv <- t(chol(vcov(fit)))
  beta.star <- beta + rv %*% rnorm(ncol(rv))
  
  # calculate psi*
  rancoef <- as.matrix(ranef(fit)[[1]]) 
  lambda <- t(rancoef)%*%rancoef
  temp <- ginv(lambda)
  ev <- eigen(temp)
  
  if (mode(ev$values) == "complex") {
    ev$values <- suppressWarnings(as.numeric(ev$values))
    ev$vectors <- suppressWarnings(matrix(as.numeric(ev$vectors),nrow=length(ev$values)))
  }
  if(sum(ev$values<0)>0)
  {
    ev$values[ev$values<0]<-0
    if(length(ev$values)>1){D<-diag(ev$values)}else{D<-ev$values}
    temp <- ev$vectors[,,drop=F]%*%D%*%t(ev$vectors[,,drop=F])
  }
  
  if(length(ev$values)>1){D<-diag(ev$values)}else{D<-ev$values}
  deco <- (ev$vectors[,,drop=F])%*%sqrt(D)
  temp.psi.star <- rWishart(1, nrow(rancoef), diag(nrow(lambda)))[,,1]
  psi.star <- ginv(deco%*%temp.psi.star%*%t(deco)) 
  
  #### psi.star positive definite?
  if (!isSymmetric(psi.star)){
    psi.star <- (psi.star + t(psi.star))/2
  }
  valprop<-eigen(psi.star)
  if(sum(valprop$values<0)>0)
  {
    valprop$values[valprop$values<0]<-0
    psi.star<-valprop$vectors%*%diag(valprop$values)%*%t(valprop$vectors)
  }
  
  
  # the main imputation task
  
  misindicator<-aggregate(as.data.frame(ry),by=list(clust=x[,clust]),FUN=function(x){sum(!x)>0})
  misindicator_sys<-aggregate(as.data.frame(ry),by=list(clust=x[,clust]),FUN=function(x){sum(x)==0})
  misindicator_spor<-aggregate(as.data.frame(ry),by=list(clust=x[,clust]),FUN=function(x){sum(x)>0})
  
  ddl<-(length(yobs)-length(fixe))#d'apres article jolani (etrange)
  sigma2<-((summary(fit)$sigma)^2*ddl)/rchisq(1,df = ddl)
  for (i in lev[misindicator[,2]]){
    # draw bi
    if(i%in%lev[misindicator_sys[,2]]){
      #systematically missing
      suppressWarnings(bi.star <- t(rmvnorm(1,mean = rep(0,nrow(psi.star)), sigma = psi.star, method="chol")))
    }
    else if(i%in%lev[misindicator_spor[,2]]){
      # sporadically missing
      Sigma2i<-sigma2*diag(sum(ry & x[,clust]==i))
      temp<-tcrossprod(psi.star,as.matrix(Z[ry & x[,clust]==i,]))%*%ginv(tcrossprod(as.matrix(Z[ry & x[,clust]==i,])%*%psi.star,as.matrix(Z[ry & x[,clust]==i,]))+Sigma2i)
      esp_bi<-temp%*%(y[ry & x[,clust]==i]-as.matrix(X[ry & x[,clust]==i,])%*%beta.star)
      var_bi<-psi.star-temp%*%as.matrix(Z[ry & x[,clust]==i,])%*%psi.star
      try(bi.star <- t(rmvnorm(1,mean = esp_bi, sigma = var_bi, method="svd")))#svd prefered to chol to avoid Nan
      if(class(bi.star)=="try-error"){
        bi.star<-esp_bi
        warning("The cov matrix for the cluster is not full rank")
      }
    }
    
    #imputation
    temp<-as.matrix(X[!ry & x[,clust]==i,]) %*% beta.star + as.matrix(Z[!ry & x[,clust]==i,])%*% bi.star
    y[!ry & x[,clust]==i]<- as.vector(temp)+rnorm(n=length(which(!ry & x[,clust]==i)),mean=0,sd=sqrt(sigma2))
  }
  return(y[!ry])
}
