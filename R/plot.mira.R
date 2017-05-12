plot.mira <-
function(x, ...){
  Mfrow_old<-par()$mfrow
  res.pool<-summary(pool(x))
  res.out<-matrix(NA,nrow = nrow(res.pool),ncol = length(x$analyses)-1)
  rownames(res.out)<-rownames(res.pool)
  colnames(res.out)<-as.character(2:length(x$analyses))
  for(ii in 2:length(x$analyses)){
    tmp<-x
    tmp$analyses<-tmp$ana[1:ii]
    res.out[,ii-1]<-2*qt(.975,summary(pool(tmp))[,"df"])*(summary(pool(tmp))[,"se"])
  }
  nbqi<-length(rownames(res.out))
  Mfrow<-c(min(nbqi,4),1+(nbqi-1)%/%4)
  par(mfrow=Mfrow)
  for(ii in seq(nrow(res.out))){
    plot(x=2:length(x$analyses),y=res.out[ii,],main=rownames(res.out)[ii],xlab="m",
                                     ylab="95% CI width")}
  par(mfrow=Mfrow_old)
}
