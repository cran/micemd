\name{find.defaultMethod}
\alias{find.defaultMethod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Suggestion of conditional imputation models to use accordingly to the incomplete dataset
}
\description{
Provides conditionnal imputation models to use for each column of the incomplete dataset according to the number of clusters, the number of individuals per cluster and the class of the variables.
}
\usage{
find.defaultMethod(don.na, ind.clust, I.small = 7, ni.small = 100, prop.small = 0.4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{don.na}{
An incomplete data frame.
}
  \item{ind.clust}{
A scalar indexes the variable corresponding to the cluster indicator.
}
  \item{I.small}{
A scalar that is used as threshold to consider the number of observed clusters (fully observed or partially observed) as small. Default is \code{I.small=7}.
}
  \item{ni.small}{
A scalar that is used as threshold to consider the number individuals per clusters (with observed values) as small. Default is \code{ni.small=100}.
}
  \item{prop.small}{
A scalar that is used as threshold to consider the number of small clusters as small. Default is \code{prop.small=0.4}.
}
}
\details{
Provides conditionnal imputation models to use for each column of the incomplete dataset according to the number of clusters, the number of individuals per cluster and the class of the variable (Audigier, V. et al 2017). Returned methods can be: 2l.stage.bin (binary), 2l.stage.norm (continuous), 2l.stage.pois (integer), 2l.glm.bin (binary), 2l.glm.norm (continuous), 2l.glm.pois (integer), 2l.jomo (continuous or binary). For a given variable, the method retained is chosen according to the following decision tree:
\tabular{rrl}{
\tab----------------------------------\tab-----------------------------------\cr
\tab Few observed \tab clusters \cr
\tab----------------------------------\tab-----------------------------------\cr
\tab Few observed values per cluster \tab Many observed values per cluster \cr
------------------\tab------------------------------------------------\tab-----------------------------------\cr
continuous	\tab 2l.glm.norm	\tab 2l.stage.norm \cr
binary	\tab 2l.glm.bin	\tab 2l.stage.bin \cr
integer	\tab 2l.glm.pois	\tab 2l.stage.pois \cr
------------------\tab------------------------------------------------\tab-----------------------------------\cr
}

\tabular{rrl}{
\tab----------------------------------\tab-----------------------------------\cr
\tab Many observed  \tab clusters\cr
\tab----------------------------------\tab-----------------------------------\cr
\tab Few observed values per cluster \tab Many observed values per cluster \cr
------------------\tab------------------------------------------------\tab-----------------------------------\cr
continuous	\tab 2l.glm.norm	\tab 2l.stage.norm\cr
binary	\tab 2l.jomo	\tab 2l.jomo\cr
integer	\tab 2l.glm.pois	\tab 2l.stage.pois\cr
------------------\tab------------------------------------------------\tab-----------------------------------\cr
}
For instance, with few observed clusters (i.e. less than \code{I.small}), and many observed values per cluster (i.e. less than \code{prop.small} clusters with less than \code{ni.small} observed values), imputation of a continuous variable according to the method 2l.stage.norm will be suggested.
}
\value{
A vector of strings with length \code{ncol(data)}.
}
\references{
Audigier, V., White, I. ,  Jolani ,S.  Debray, T.,  Quartagno, M., Carpenter, J., van Buuren, S. and Resche-Rigon, M.
Multiple imputation for multilevel data with continuous and binary variables. ArXiv e-prints <arXiv:1702.00971>, 2017.


Jolani, S., Debray, T. P. A., Koffijberg, H., van Buuren, S., and Moons, K. G. M. (2015). Imputation
of systematically missing predictors in an individual participant data meta-analysis: a
generalized approach using MICE. Statistics in Medicine, 34(11):1841{1863}. <doi:10.1002/sim.6451>


Quartagno, M. and Carpenter, J. R. (2016). Multiple imputation for IPD meta-analysis: allowing
for heterogeneity and studies with missing covariates. Statistics in Medicine, 35(17):2938{
2954}. <doi:10.1002/sim.6837>

Resche-Rigon, M. and White, I. R. (2016). Multiple imputation by chained equations for systematically
and sporadically missing multilevel data. Statistical Methods in Medical Research.
To appear. <doi:10.1177/0962280216666564>
}
\author{
Vincent Audigier \email{vincent.audigier@univ-paris-diderot.fr}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mice}}, \code{\link{mice.par}}
}
\examples{

  data(CHEM97Na)
  
  ind.clust<-1#index for the cluster variable
  
  #initialisation of the argument predictorMatrix
  temp<-mice(CHEM97Na,m=1,maxit=0)
  temp$pred[ind.clust,ind.clust]<-0;temp$pred[-ind.clust,ind.clust]<- -2
  temp$pred[temp$pred==1]<-2
  predictor.matrix<-temp$pred
  
  #initialisation of the argument method
  method<-find.defaultMethod(CHEM97Na,ind.clust)
  print(method)
  
  #multiple imputation by chained equations (parallel calculation)
  #res.mice<-mice.par(CHEM97Na,m=3,predictorMatrix = predictor.matrix,method=method)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{mice}
