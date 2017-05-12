mice.par <-
function(don.na, m = 5, method = vector("character", length = ncol(don.na)), 
                   predictorMatrix = (1 - diag(1, ncol(don.na))),
                   visitSequence = (1:ncol(don.na))[apply(is.na(don.na),2, any)], form = vector("character", length = ncol(don.na)), 
                   post = vector("character", length = ncol(don.na)), defaultMethod = c("pmm", 
                                                                                        "logreg", "polyreg", "polr"), maxit = 5, diagnostics = TRUE, 
                   seed = NA, imputationMethod = NULL, defaultImputationMethod = NULL, 
                   data.init = NULL, nnodes=5, path.outfile=NULL,...){

  cl <- makeCluster(nnodes, type="PSOCK")


  tmp<-list(...)
  if("k"%in%names(tmp)){k<-tmp[["k"]]}else{k<-5}
  if("method_est"%in%names(tmp)){method_est<-tmp[["method_est"]]}else{method_est<-"mm"}
  if("incluster"%in%names(tmp)){incluster<-tmp[["incluster"]]}else{incluster<-FALSE}
  if("nburn"%in%names(tmp)){nburn<-tmp[["nburn"]]}else{nburn<-200}

  

  clusterExport(cl, list("mice.impute.2l.glm.bin",
                         "mice.impute.2l.glm.norm",
                         "mice.impute.2l.2stage.bin",
                         "mice.impute.2l.2stage.bin.intern",
                         "mice.impute.2l.2stage.norm.intern",
                         "mice.impute.2l.2stage.norm",
                         "mice.impute.2l.2stage.pmm",
                         "mice.impute.2l.jomo",
                         "mice.impute.2l.2stage.pois.intern",
                         "mice.impute.2l.2stage.pois",
                         "mice.impute.2l.glm.pois",
                         "mice.impute.2l.norm", "mice.impute.2l.pan", "mice.impute.2lonly.mean", 
                         "mice.impute.2lonly.norm", "mice.impute.2lonly.pmm", "mice.impute.cart", 
                         "mice.impute.fastpmm", "mice.impute.lda", "mice.impute.logreg", 
                         "mice.impute.logreg.boot", "mice.impute.mean", "mice.impute.midastouch", 
                         "mice.impute.norm", "mice.impute.norm.boot", "mice.impute.norm.nob", 
                         "mice.impute.norm.predict", "mice.impute.passive", "mice.impute.pmm", 
                         "mice.impute.polr", "mice.impute.polyreg", "mice.impute.quadratic", 
                         "mice.impute.rf", "mice.impute.ri", "mice.impute.sample",
                         "don.na","method",
                         "predictorMatrix", 
                         "visitSequence",
                         "post",
                         "defaultMethod",
                         "maxit",
                         "diagnostics", 
                         "imputationMethod",
                         "defaultImputationMethod", 
                         "data.init",
                         "maxit",
                         "find.defaultMethod",
                         "nnodes",
                         "k",
                         "method_est",
                         "incluster",
                         "nburn","path.outfile"),envir = environment())  
  if(!is.null(path.outfile)){
    clusterEvalQ(cl, sink(paste0(path.outfile,"/output", Sys.getpid(), ".txt")))
  }
  
  
  res<-parSapply(cl,as.list(1:m), FUN=function(mtmp,don.na,method,predictorMatrix,visitSequence,
                                               post,imputationMethod,defaultImputationMethod,data.init,
                                               maxit,m,nnodes,k,method_est,incluster,nburn){
    hexval <- paste0("0x", sapply(1:m, digest, "crc32")) 
    seedlist <- type.convert(hexval) %% .Machine$integer.max
    res.mice<-mice(data=don.na,m=1,method = method,predictorMatrix=predictorMatrix,
                   visitSequence=visitSequence,post=post,imputationMethod=imputationMethod,
                   defaultImputationMethod=defaultImputationMethod,data.init=data.init,
                   maxit = maxit,printFlag=TRUE,seed=seedlist[mtmp],k=k,method_est=method_est,incluster=incluster,nburn=nburn)
  },don.na=don.na,method=method,predictorMatrix=predictorMatrix,visitSequence=visitSequence,post=post,
  imputationMethod=imputationMethod,defaultImputationMethod=defaultImputationMethod,data.init=data.init,
  maxit=maxit,m=m,nnodes=nnodes,k=k,method_est=method_est,incluster=incluster,nburn=nburn,simplify = FALSE)

  stopCluster(cl)
  
  res.out<-res[[1]]
  res.out$call<-  match.call()
  res.out$m<-m
  res.out$imp<-mapply(as.list(colnames(don.na)),FUN=function(xx,res){do.call(cbind,lapply(lapply(res,"[[","imp"),"[[",xx))},MoreArgs=list(res=res))
  names(res.out$imp)<-colnames(don.na)
  res.out$imp<-lapply(res.out$imp,function(xx){if(!is.null(xx)){yy<-xx;colnames(yy)<-as.character(seq(ncol(xx)));return(yy)}else{return(xx)}})
  res.out$seed<-unlist(lapply(res,"[[","seed"))
  res.out$lastSeedValue<-lapply(res,"[[","lastSeedValue")
  res.out$chainMean<-do.call(abind,lapply(res,"[[","chainMean"),3)
  dimnames(res.out$chainMean)[[3]]<-paste("Chain",seq(m))
  res.out$chainVar<-do.call(abind,lapply(res,"[[","chainVar"),3)
  dimnames(res.out$chainVar)[[3]]<-paste("Chain",seq(m))
  res.out$loggedEvents<-lapply(res,"[[","loggedEvents")
  class(res.out)<-"mids"
  return(res.out)
}
