overimpute <-
function (res.mice, plotvars = NULL,plotinds=NULL,nnodes=5,path.outfile =NULL,alpha=0.1) 
{
  if (res.mice$m < 100) {
    warning("The number of imputed data sets is too low to build confidence intervals according to the quantiles method. You should run mice with m over than 100.")
  }
 
  if (is.null(plotvars)) {
    Var <- seq(ncol(res.mice$data))
  }
  else {
    Var <- plotvars
  }
  if (is.null(plotinds)) {
    Ind <- seq(nrow(res.mice$data))
  }
  else {
    Ind <- plotinds
  }
  
  #index for cluster
  tmp<-apply(res.mice$predictorMatrix,2,FUN=function(xx){any(-2==xx)})
  if(any(tmp)){
  ind.clust<-which(names(which(tmp))%in%colnames(res.mice$data))
  }else{
    ind.clust<-0
  }
  
  #no plot for categorical variables
  is.plot<-sapply(res.mice$data[,Var],is.numeric)
  #no plot for binary variable
  is.plot[apply(res.mice$data[,Var],2,function(xx){length(table(xx))==2})]<-FALSE
  #no plot for cluster index
  is.plot[ind.clust]<-FALSE
  #no plot for variables without any missing values
  is.plot[colSums(is.na(res.mice$data[,Var]))==0]<-FALSE
  
  don.plot<-as.matrix(res.mice$data[Ind,is.plot])
  res.over<-matrix(NA, nrow=sum(!is.na(don.plot)),ncol=res.mice$m)
  
   cl <- makeCluster(nnodes,type = "PSOCK")
   clusterExport(cl, list("is.plot","res.mice","path.outfile","Ind"), envir = environment())
   if (!is.null(path.outfile)) {
     clusterEvalQ(cl, sink(paste0(path.outfile, "/output", 
                                  Sys.getpid(), ".txt")))
   }
  clusterEvalQ(cl, library(mice))
  clusterEvalQ(cl, library(micemd))
  
  res.over<-t(parSapply(cl,which(!is.na(don.plot)),FUN=function(jj,is.plot,res.mice,Ind){
    cat("cell number : ",jj,"\n")
    sapply(seq(res.mice$m),FUN=function(m,ii,is.plot,res.mice,Ind){
      don.over<-complete(res.mice,m)
      tmp<-as.matrix(don.over[Ind,names(which(is.plot))])
      tmp[ii]<-NA
      don.over[Ind,names(which(is.plot))]<-tmp
      res.mice.tmp<-try(mice(don.over,m=1,maxit=1,printFlag = FALSE,method = res.mice$method, predictorMatrix = res.mice$predictorMatrix))
      if(class(res.mice.tmp)=="try-error"){
        res<-NA
        }else{
      res<-unlist(res.mice.tmp$imp)
      }
      return(res)
    },ii=jj,is.plot=is.plot,res.mice=res.mice,Ind=Ind)},is.plot=is.plot,res.mice=res.mice,Ind=Ind))
  stopCluster(cl)
  
  res.plot <- t(apply(res.over, 1, function(x) {
    xbar <- mean(x)
    temp <- quantile(x, probs = c(alpha/2, 1-alpha/2),na.rm=TRUE)
    binf <- temp[[1]]
    bsup <- temp[[2]]
    return(c(xbar = xbar, binf = binf, bsup = bsup))
  }))
  
  pct <- rep(rowMeans(is.na(res.mice$data[Ind,])),ncol(don.plot))
  col <- cut(pct, c(-0.1, 0.2, 0.4, 0.6, 0.8, 1.1))
  levels(col) <- c("blue", "green", heat.colors(3)[c(3, 
                                                     2, 1)])
  col <- as.character(col)
  res.over.value<-res.over
  res.over<-cbind.data.frame(var=unlist(mapply(FUN = rep,names(is.plot[is.plot]),each=apply(!is.na(don.plot),2,sum))),
                              trueval=don.plot[!is.na(don.plot)],res.plot,pct=pct[!is.na(don.plot)],col=col[!is.na(don.plot)])
  
  par(mfrow = c(ceiling(sqrt(length(is.plot[is.plot]))), ceiling(length(is.plot[is.plot])/ceiling(sqrt(length(is.plot[is.plot]))))), 
      mar = c(5, 4, 4, 2) - 1.9)
  by(res.over,INDICES =  res.over$var,FUN=function(xx){
    plot(x = xx[,"trueval"], y = xx[, "xbar"], col = as.character(xx[,"col"]), 
         xlab = "observed values", ylab = "imputed values", 
         main = xx[1,"var"], ylim = c(min(xx[, 
                                             "binf"], na.rm = T), max(xx[, "bsup"], na.rm = T)))
    abline(0, 1)
    segments(x0 = xx[,"trueval"], x1 = xx[,"trueval"], 
             y0 = xx[, "binf"], y1 = xx[, "bsup"], col =  as.character(xx[,"col"]))
    legend("topleft", legend = c("0-0.2", "0.2-0.4", "0.4-0.6", 
                                 "0.6-0.8", "0.8-1"), col = c("blue", "green", heat.colors(3)[c(3, 
                                                                                                2, 1)]), bty = "n", lty = 1, horiz = F, cex = 0.7, 
           lwd = 0.4)
  }
  )
  return(list(res.plot=res.over,res.values=res.over.value))
}
