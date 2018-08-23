mice.impute.2l.jomo <-
function(y, ry, x,type,nburn=200,...){
  finddata<-function(y,ry,x,type){
    
    #identification indicatrice
    debut<-which(regexpr(pattern = ".1",text = names(x))>0)
    varquali<-gsub(".1",replacement = "",names(x)[debut])
    indexindic<-sapply(varquali,FUN=function(nom,xx){which(regexpr(pattern = nom,text = names(xx))>0)},xx=x)
    nomvar<-names(x);nomvar[unlist(indexindic)]<-NA;nomvar[is.na(nomvar)]<-names(indexindic)
    names(type)<-unique(nomvar)
    #transformation en facteur
    res<-sapply(indexindic,FUN=function(variable,xx){
      tab.0<-1-rowSums(xx[,variable,drop=FALSE])
      tab<-cbind.data.frame(tab.0,xx[,variable])
      res<-apply(tab,1,which.max)
      res<-sapply(res,as.factor)
      return(res)
    },xx=x)
    if(length(indexindic)>0){
      don<-cbind.data.frame(x[,-unlist(indexindic),drop=FALSE],res)
    }else{
      don<-x
    }
    res<-list(x=don,type=type[names(don)])#clust,quanti, quali
    return(res)
  }
  
  jomo1ran.silent<-function (Y, X = NULL, Z = NULL, clus, beta.start = NULL, u.start = NULL, 
                             l1cov.start = NULL, l2cov.start = NULL, l1cov.prior = NULL, 
                             l2cov.prior = NULL, nburn = nburn, nbetween = 1, nimp = 2, 
                             a = NULL, a.prior = NULL, meth = "random", output = 0, out.iter = 1000) {
    
    stopifnot(meth == "common" | meth == "fixed" | meth == "random")
    ncon = 0
    ncat = 0
    Y.con = NULL
    Y.cat = NULL
    Y.numcat = NULL
    for (i in 1:ncol(Y)) {
      if (is.numeric(Y[, i])) {
        ncon = ncon + 1
        if (is.null(Y.con)) {
          Y.con <- data.frame(Y[, i])
        }
        else {
          Y.con <- data.frame(Y.con, Y[, i])
        }
        colnames(Y.con)[ncon] <- colnames(Y)[i]
      }
      else {
        if (is.factor(Y[, i])) {
          ncat = ncat + 1
          if (is.null(Y.cat)) {
            Y.cat <- data.frame(Y[, i])
          }
          else {
            Y.cat <- data.frame(Y.cat, Y[, i])
          }
          colnames(Y.cat)[ncat] <- colnames(Y)[i]
          Y.numcat <- cbind(Y.numcat, nlevels(Y[, i]))
        }
      }
    }
    if (is.null(X)) 
      X = matrix(1, nrow(Y), 1)
    if (is.null(Z)) 
      Z = matrix(1, nrow(Y), 1)
    if (meth == "common") {
      if (ncat == 0 & ncon > 0) {
        # cat("Found ", ncon, "continuous outcomes and no categorical. Using function jomo1rancon.", 
        # "\n")
        imp <- jomo1rancon(Y = Y.con, X = X, Z = Z, clus = clus, 
                           beta.start = beta.start, u.start = u.start, l1cov.start = l1cov.start, 
                           l2cov.start = l2cov.start, l1cov.prior = l1cov.prior, 
                           l2cov.prior = l2cov.prior, nburn = nburn, nbetween = nbetween, 
                           nimp = nimp, output = output, out.iter = out.iter)
      }
      if (ncat > 0 & ncon == 0) {
        # cat("Found ", ncat, "categorical outcomes and no continuous. Using function jomo1rancat.", 
        # "\n")
        imp <- jomo1rancat(Y.cat = Y.cat, Y.numcat = Y.numcat, 
                           X = X, Z = Z, clus = clus, beta.start = beta.start, 
                           u.start = u.start, l1cov.start = l1cov.start, 
                           l2cov.start = l2cov.start, l1cov.prior = l1cov.prior, 
                           l2cov.prior = l2cov.prior, nburn = nburn, nbetween = nbetween, 
                           nimp = nimp, output = output, out.iter = out.iter)
      }
      if (ncat > 0 & ncon > 0) {
        #         cat("Found ", ncon, "continuous outcomes and ", ncat, 
        #             "categorical. Using function jomo1ranmix.", "\n")
        imp <- jomo1ranmix(Y.con, Y.cat = Y.cat, Y.numcat = Y.numcat, 
                           X = X, Z = Z, clus = clus, beta.start = beta.start, 
                           u.start = u.start, l1cov.start = l1cov.start, 
                           l2cov.start = l2cov.start, l1cov.prior = l1cov.prior, 
                           l2cov.prior = l2cov.prior, nburn = nburn, nbetween = nbetween, 
                           nimp = nimp, output = output, out.iter = out.iter)
      }
    }
    if (meth == "fixed") {
      if (ncat == 0 & ncon > 0) {
        #         cat("Found ", ncon, "continuous outcomes and no categorical. Using function jomo1ranconhr with fixed cluster-specific covariance matrices.", 
        #             "\n")
        imp <- jomo1ranconhr(Y = Y.con, X = X, Z = Z, clus = clus, 
                             beta.start = beta.start, u.start = u.start, l1cov.start = l1cov.start, 
                             l2cov.start = l2cov.start, l1cov.prior = l1cov.prior, 
                             l2cov.prior = l2cov.prior, nburn = nburn, nbetween = nbetween, 
                             nimp = nimp, a = 0, meth = "fixed", output = output, 
                             out.iter = out.iter)
      }
      if (ncat > 0 & ncon == 0) {
        #         cat("Found ", ncat, "categorical outcomes and no continuous. Using function jomo1rancathr with fixed cluster-specific covariance matrices.", 
        #             "\n")
        imp <- jomo1rancathr(Y.cat = Y.cat, Y.numcat = Y.numcat, 
                             X = X, Z = Z, clus = clus, beta.start = beta.start, 
                             u.start = u.start, l1cov.start = l1cov.start, 
                             l2cov.start = l2cov.start, l1cov.prior = l1cov.prior, 
                             l2cov.prior = l2cov.prior, nburn = nburn, nbetween = nbetween, 
                             nimp = nimp, a = 0, meth = "fixed", output = output, 
                             out.iter = out.iter)
      }
      if (ncat > 0 & ncon > 0) {
        #         cat("Found ", ncon, "continuous outcomes and ", ncat, 
        #             "categorical. Using function jomo1ranmixhr with fixed cluster-specific covariance matrices.", 
        #             "\n")
        imp <- jomo1ranmixhr(Y.con, Y.cat = Y.cat, Y.numcat = Y.numcat, 
                             X = X, Z = Z, clus = clus, beta.start = beta.start, 
                             u.start = u.start, l1cov.start = l1cov.start, 
                             l2cov.start = l2cov.start, l1cov.prior = l1cov.prior, 
                             l2cov.prior = l2cov.prior, nburn = nburn, nbetween = nbetween, 
                             nimp = nimp, a = 0, meth = "fixed", output = output, 
                             out.iter = out.iter)
      }
    }
    if (meth == "random") {
      if (ncat == 0 & ncon > 0) {
        #         cat("Found ", ncon, "continuous outcomes and no categorical. Using function jomo1ranconhr with random cluster-specific covariance matrices.", 
        #             "\n")
        imp <- jomo1ranconhr(Y = Y.con, X = X, Z = Z, clus = clus, 
                             beta.start = beta.start, u.start = u.start, l1cov.start = l1cov.start, 
                             l2cov.start = l2cov.start, l1cov.prior = l1cov.prior, 
                             l2cov.prior = l2cov.prior, nburn = nburn, nbetween = nbetween, 
                             nimp = nimp, a = a, a.prior = a.prior, meth = "random", 
                             output = output, out.iter = out.iter)
      }
      if (ncat > 0 & ncon == 0) {
        #         cat("Found ", ncat, "categorical outcomes and no continuous. Using function jomo1rancathr with random cluster-specific covariance matrices.", 
        #             "\n")
        imp <- jomo1rancathr(Y.cat = Y.cat, Y.numcat = Y.numcat, 
                             X = X, Z = Z, clus = clus, beta.start = beta.start, 
                             u.start = u.start, l1cov.start = l1cov.start, 
                             l2cov.start = l2cov.start, l1cov.prior = l1cov.prior, 
                             l2cov.prior = l2cov.prior, nburn = nburn, nbetween = nbetween, 
                             nimp = nimp, a = a, a.prior = a.prior, meth = "random", 
                             output = output, out.iter = out.iter)
      }
      if (ncat > 0 & ncon > 0) {
        #         cat("Found ", ncon, "continuous outcomes and ", ncat, 
        #             "categorical. Using function jomo1ranmixhr with random cluster-specific covariance matrices.", 
        #             "\n")
        imp <- jomo1ranmixhr(Y.con, Y.cat = Y.cat, Y.numcat = Y.numcat, 
                             X = X, Z = Z, clus = clus, beta.start = beta.start, 
                             u.start = u.start, l1cov.start = l1cov.start, 
                             l2cov.start = l2cov.start, l1cov.prior = l1cov.prior, 
                             l2cov.prior = l2cov.prior, nburn = nburn, nbetween = nbetween, 
                             nimp = nimp, a = a, a.prior = a.prior, meth = "random", 
                             output = output, out.iter = out.iter)
      }
    }
    return(imp)
    
  }
  
  xx<-as.data.frame(finddata(y,ry,x,type)$x)
  names(xx) <- paste("V",1:ncol(x),sep="")
  clust <- "V1"
  fixe <- names(xx)[type>0]
  yna<-y
  yna[!ry]<-NA
  
  res.jomo<-try(jomo1ran.silent(Y=cbind.data.frame(yna,xx[,fixe,drop=FALSE]),clus=xx[,clust,drop=FALSE],
                                nburn=nburn,
                                nbetween=1, nimp=2, output=0, out.iter=nburn,meth="random"))
  
  if(class(res.jomo)!="try-error"){
    ind.y<-apply(res.jomo,2,FUN=function(xx){sum(is.na(xx))>0})
    res.out<-res.jomo[which(res.jomo[,"Imputation"]==1),ind.y][!ry]
    return(res.out)
  }else{warning("jomo, imputation fails");return(y[!ry])}
}
