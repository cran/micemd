find.defaultMethod <-
function(don.na,ind.clust,I.small=7,ni.small=100,prop.small=0.4){
  #don.na : an incomplete dataset with the cluster variable in first column

  tmp<-mapply(as.data.frame(don.na),FUN=function(xx,v.clust){
    piNA<-as.vector(by(xx,INDICES = v.clust,FUN = function(yy){sum(is.na(yy))/length(yy)}))#proportion of NA per cluster
    nbobs<-as.vector(by(xx,INDICES = v.clust,FUN = function(yy){sum(!is.na(yy))}))#number of observed values per cluster

    if(sum(is.na(xx))==0){return("")}
    
    if(sum(piNA!=0)<=I.small){
      #the number of observed clusters is small
      if(mean(nbobs[nbobs>0]<=ni.small)>prop.small){
        #the number of observed values per cluster is often small
        if(is.factor(xx)){return("2l.glm.bin")}else if(is.integer(xx)){return("2l.glm.pois")}else{return("2l.glm.norm")}
      }else{
        #the number of observed values per cluster is generally not small
        if(is.factor(xx)){return("2l.2stage.bin")}else if(is.integer(xx)){return("2l.2stage.pois")}else{return("2l.2stage.norm")}
      }
      
    }else{
      #the number of observed clusters large
      if(mean(nbobs[nbobs>0]<=ni.small)>prop.small){
        #the number of observed values per cluster is often small
        if(is.factor(xx)){return("2l.jomo")}else if(is.integer(xx)){return("2l.glm.pois")}else{return("2l.glm.norm")}
      }else{
        #the number of observed values per cluster is generally not small
        if(is.factor(xx)){return("2l.jomo")}else if(is.integer(xx)){return("2l.2stage.pois")}else{return("2l.2stage.norm")}
      }
    }
  },MoreArgs = list(v.clust=don.na[,ind.clust]))
  return(tmp)
}
