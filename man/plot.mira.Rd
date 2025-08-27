\name{plot.mira}
\alias{plot.mira}
\title{
Graphical investigation for the number of generated datasets
}
\description{
The \code{plot} method for a \code{mira} object plots the confidence interval length against the number of multiply imputed datasets from 2 to \code{m}. This is a graphical tool to check if the variability due to the simulation of the multiple imputation process can be substantially reduced by increasing the number of generated datasets \code{m}.
}
\usage{
\method{plot}{mira}(x, \dots)
}

\arguments{
  \item{x}{
An object of class \code{mira}.
}
  \item{\dots}{
Extra arguments for \code{\link{plot.mira}}
}
}


\references{
Schafer, J. L. (1997). Analysis of Incomplete Multivariate Data. Chapman & Hall/CRC, London
}
\author{
Vincent Audigier \email{vincent.audigier@cnam.fr}
}

\seealso{
\code{\link[mice]{mice}} \code{\link[mice]{mids}}
}
\examples{

  require(nlme)
  data(CHEM97Na)
  
  ind.clust<-1#index for the cluster variable
  
  #initialisation of the argument predictorMatrix
  predictor.matrix<-mice(CHEM97Na,m=1,maxit=0)$pred
  predictor.matrix[ind.clust,ind.clust]<-0
  predictor.matrix[-ind.clust,ind.clust]<- -2
  predictor.matrix[predictor.matrix==1]<-2

  
  #initialisation of the argument method
  method<-c("", "2l.2stage.bin", "2l.2stage.pois", "2l.2stage.norm", "") #quickest methods
  
  #multiple imputation by chained equations (parallel calculation)
  #res.mice<-mice.par(CHEM97Na,m=15,predictorMatrix = predictor.matrix,method=method)

  #analysis (apply a linear mixed effects model to each imputed dataset)
  #ana<-with(res.mice,expr=lme(fixed=formula(Score~Sex+GSCE+Age),
  #                            random=formula(~1|School),method="REML",
  #                            control=list(maxIter=100,msMaxIter=100,niterEM=25)))

  #graphical investigation for the number of generated datasets m
  #plot(ana)

}

