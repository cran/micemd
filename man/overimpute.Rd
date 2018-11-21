\name{overimpute}
\alias{overimpute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Overimputation diagnostic plot
}
\description{
Assess the fit of the predictive distribution after performing multiple imputation with mice
}
\usage{
overimpute(res.mice, plotvars = NULL, plotinds = NULL,
  nnodes = 5, path.outfile = NULL, alpha = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res.mice}{
An object of class mids
}
  \item{plotvars}{
column index of the variables overimputed
}
  \item{plotinds}{
row index of the individuals overimputed

}
  \item{nnodes}{
A scalar indicating the number of nodes for parallel calculation. Default value is 5.
}
  \item{path.outfile}{
A vector of strings indicating the path for redirection of print messages. Default value is NULL, meaning that silent imputation is performed. Otherwise, print messages are saved in the files path.outfile/output.txt. One file per node is generated.
}
  \item{alpha}{
alpha level for prediction intervals
}
}
\details{
This function imputes each observed values from each of the parameters of the imputation model obtained from the mice procedure. The comparison between the "overimputed" values and the observed values is made by building a confidence interval for each observed value using the quantiles of the overimputed values (Blackwell et al. (2015)). Note that confidence intervals builded with quantiles require a large number of imputations. If the model fits well the data, then the 90\% confidence interval should contain the observed value in 90\% of the cases. The function Overimpute takes as an input the output of the mice or mice.par function (res.mice), the indices of the incomplete continuous variables that are plotted (plotvars), the indices of individuals (can be useful for time consumming imputation methods), the number of nodes for parallel computation, and the path for exporting print message generated during the parallel process.
}
\value{
A list of two matrices
 \item{res.plot}{7-columns matrix that contains (1) the variable which is overimputed, (2) the observed value of the observation, (3) the mean of the overimputations, (4) the lower bound of the confidence interval of the overimputations, (5) the upper bound of the confidence interval of the overimputations, (6) the proportion of the other variables that were missing for that observation in the original data, and (7) the color for graphical representation.}
 \item{res.values}{A matrix with overimputed values for each cell. The number of columns corresponds to the number of values generated (i.e. the number of imputed tables)}

}
\references{
Blackwell, M., Honaker, J. and King. G. 2015. A Unified Approach to Measurement Error and Missing Data: Overview and Applications. Sociological Methods and Research, 1-39. <doi:10.1177/0049124115585360>
}
\author{
Vincent Audigier \email{vincent.audigier@cnam.fr}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mice},\link{parallel}, \link{mice.par}}
}

\examples{

require(parallel)
nnodes<-detectCores()-1#number of nodes
m<-100#nb generated values per observation

################
#one level data
################
require(mice)
data(nhanes)
#res.mice<-mice.par(nhanes,m = m,nnodes = nnodes)
#res.over<-overimpute(res.mice, nnodes = nnodes)

################
#two level data (time consumming)
################
data(CHEM97Na)

ind.clust<-1#index for the cluster variable

#initialisation of the argument predictorMatrix
temp<-mice(CHEM97Na,m=1,maxit=0)
temp$pred[ind.clust,ind.clust]<-0
temp$pred[-ind.clust,ind.clust]<- -2
predictor.matrix<-temp$pred

#initialisation of the argument method
method<-find.defaultMethod(CHEM97Na,ind.clust)

#multiple imputation by chained equations (time consumming)

#res.mice<-mice.par(CHEM97Na,
#                   predictorMatrix = predictor.matrix,
#                   method=method,m=m,nnodes = nnodes)


#overimputation on 30 individuals
#res.over<-overimpute(res.mice,
#                     nnodes=nnodes,
#                     plotinds=sample(x = seq(nrow(CHEM97Na)),size = 30))

}
\keyword{mice}
