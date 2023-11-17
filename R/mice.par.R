mice.par <-
function(don.na, m = 5, method = NULL, predictorMatrix= (1 - diag(1, ncol(don.na))), where = NULL, 
                   visitSequence = NULL, blots = NULL, post = NULL, blocks,formulas,
                   defaultMethod = c("pmm", "logreg", "polyreg", "polr"), maxit = 5, 
                  seed = NA, data.init = NULL, nnodes=5, path.outfile=NULL,...){
  
  
  
  cl <- makeCluster(nnodes, type="PSOCK")
  if(!is.na(seed)){clusterSetRNGStream(cl,seed)}
  
  tmp<-list(...)
  if("kpmm"%in%names(tmp)){kpmm<-tmp[["kpmm"]]}else{kpmm<-5}
  if("method_est"%in%names(tmp)){method_est<-tmp[["method_est"]]}else{method_est<-"mm"}
  if("incluster"%in%names(tmp)){incluster<-tmp[["incluster"]]}else{incluster<-FALSE}
  if("nburn"%in%names(tmp)){nburn<-tmp[["nburn"]]}else{nburn<-200}
  if("ypmm"%in%names(tmp)){ypmm<-tmp[["ypmm"]]}else{ypmm<-NULL}
  if("pmm"%in%names(tmp)){pmm<-tmp[["pmm"]]}else{pmm<-FALSE}
  if("pred_std"%in%names(tmp)){pred_std<-tmp[["pred_std"]]}else{pred_std<-FALSE}
  if("meta_method"%in%names(tmp)){meta_method<-tmp[["meta_method"]]}else{meta_method <- "reml"}
  
  
  if(maxit==0){stop("The argument maxit=0 is not relevant for parallel calculation, use the mice function from the mice package")}
  
  
  clusterExport(cl, list("mice","mice.impute.2l.glm.bin","mice.impute.2l.2stage.heckman",
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
                         "mice.impute.2l.lmer",
                         "mice.impute.2l.norm", "mice.impute.2l.pan", "mice.impute.2lonly.mean", 
                         "mice.impute.2lonly.norm", "mice.impute.2lonly.pmm", "mice.impute.cart", 
                         "mice.impute.lda", "mice.impute.logreg", 
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
                         "where", 
                         "blots", 
                         "data.init",
                         "maxit",
                         "find.defaultMethod",
                         "nnodes",
                         "kpmm",
                         "method_est",
                         "incluster",
                         "nburn","path.outfile","ypmm","pmm","pred_std","meta_method"),envir = environment())  
  
  if(!missing(blocks)){clusterExport(cl, list("blocks"),envir = environment())}
  if(!missing(formulas)){clusterExport(cl, list("formulas"),envir = environment())} 
  
  if(!is.null(path.outfile)){
    clusterEvalQ(cl, sink(paste0(path.outfile,"/output", Sys.getpid(), ".txt")))
  }
  
  
  if(!missing(blocks)&!missing(formulas)){
    res<-parSapply(cl,as.list(1:m), FUN=function(mtmp,don.na,method,predictorMatrix,visitSequence,where,blocks,formulas,blots,
                                                 post,defaultMethod,data.init,
                                                 maxit,nnodes,kpmm,method_est,incluster,nburn,ypmm,pmm,pred_std,meta_method){
      
      res.mice<-mice(data=don.na,m=1,method = method,predictorMatrix=predictorMatrix,
                     visitSequence=visitSequence,where=where,blocks=blocks,formulas=formulas,blots=blots,post=post,
                     defaultMethod=defaultMethod,data.init=data.init,
                     maxit = maxit,printFlag=TRUE,seed=NA,method_est=method_est,incluster=incluster,nburn=nburn,kpmm=kpmm,ypmm=ypmm,pmm=pmm,pred_std=pred_std,meta_method=meta_method)
    },don.na=don.na,
    method=method,
    predictorMatrix=predictorMatrix,
    visitSequence=visitSequence,
    blocks=blocks,
    formulas=formulas,
    where=where,
    blots=blots,
    post=post,
    defaultMethod=defaultMethod,
    data.init=data.init,
    maxit=maxit,
    nnodes=nnodes,
    method_est=method_est,
    kpmm=kpmm,
    ypmm=ypmm,
    pmm=pmm,
    pred_std=pred_std,
    meta_method=meta_method,
    incluster=incluster,
    nburn=nburn,
    simplify = FALSE)
  }else if(missing(blocks)&missing(formulas)){
    res<-parSapply(cl,as.list(1:m), FUN=function(mtmp,don.na,method,predictorMatrix,visitSequence,where,blots,
                                                 post,defaultMethod,data.init,
                                                 maxit,nnodes,kpmm,method_est,incluster,nburn,ypmm,pmm,pred_std,meta_method){
      
      res.mice<-mice(data=don.na,m=1,method = method,predictorMatrix=predictorMatrix,
                     visitSequence=visitSequence,where=where,blots=blots,post=post,
                     defaultMethod=defaultMethod,data.init=data.init,
                     maxit = maxit,printFlag=TRUE,seed=NA,method_est=method_est,incluster=incluster,nburn=nburn,kpmm=kpmm,ypmm=ypmm,pmm=pmm,pred_std=pred_std,meta_method=meta_method)
    },don.na=don.na,
    method=method,
    predictorMatrix=predictorMatrix,
    visitSequence=visitSequence,
    where=where,
    blots=blots,
    post=post,
    defaultMethod=defaultMethod,
    data.init=data.init,
    maxit=maxit,
    nnodes=nnodes,
    method_est=method_est,
    kpmm=kpmm,
    ypmm=ypmm,
    pmm=pmm,
    pred_std=pred_std,
    meta_method=meta_method,
    incluster=incluster,
    nburn=nburn,
    simplify = FALSE)
    
  }else{stop("blocks or formulas arguments are not defined. Currently, this case is not handled by the mice.par function")}
  
  stopCluster(cl)
  
  res.out<-res[[1]]
  res.out$call<-  match.call()
  res.out$m<-m
  res.out$imp<-mapply(as.list(colnames(don.na)),FUN=function(xx,res){do.call(cbind,lapply(lapply(res,"[[","imp"),"[[",xx))},MoreArgs=list(res=res),SIMPLIFY = FALSE)
  names(res.out$imp)<-colnames(don.na)
  res.out$imp<-lapply(res.out$imp,function(xx){if(!is.null(xx)){yy<-xx;colnames(yy)<-as.character(seq(ncol(xx)));return(yy)}else{return(xx)}})
  res.out$seed<-seed
  res.out$lastSeedValue<-lapply(res,"[[","lastSeedValue")
  res.out$chainMean<-do.call(abind,lapply(res,"[[","chainMean"),3)
  dimnames(res.out$chainMean)[[3]]<-paste("Chain",seq(m))
  res.out$chainVar<-do.call(abind,lapply(res,"[[","chainVar"),3)
  dimnames(res.out$chainVar)[[3]]<-paste("Chain",seq(m))
  res.out$loggedEvents<-lapply(res,"[[","loggedEvents")
  class(res.out)<-"mids"
  return(res.out)
}
