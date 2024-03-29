\name{mice.par}
\alias{mice.par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parallel calculations for Multivariate Imputation by Chained Equations
}
\description{
Parallel calculations for Multivariate Imputation by Chained Equations using the R package \code{parallel}.
}
\usage{
mice.par(don.na, m = 5, method = NULL, predictorMatrix, where = NULL,
visitSequence = NULL, blots = NULL, post = NULL, blocks, formulas,
defaultMethod = c("pmm", "logreg", "polyreg", "polr"), maxit = 5,
seed = NA, data.init = NULL, nnodes = 5, path.outfile = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{don.na}{
A data frame or a matrix containing the incomplete data.  Missing
values are coded as \code{NA}.
}
  \item{m}{
Number of multiple imputations. The default is \code{m=5}.
}
  \item{method}{
Can be either a single string, or a vector of strings with
length \code{ncol(data)}, specifying the elementary imputation method to be
used for each column in data. If specified as a single string, the same
method will be used for all columns.  The default imputation method (when no
argument is specified) depends on the measurement level of the target column
and are specified by the \code{defaultMethod} argument.  Columns that need
not be imputed have the empty method \code{''}.  See details for more
information.
}
  \item{predictorMatrix}{
A square matrix of size \code{ncol(data)} containing
0/1 data specifying the set of predictors to be used for each target column.
Rows correspond to target variables (i.e. variables to be imputed), in the
sequence as they appear in data. A value of '1' means that the column
variable is used as a predictor for the target variable (in the rows). The
diagonal of \code{predictorMatrix} must be zero. The default for
\code{predictorMatrix} is that all other columns are used as predictors
(sometimes called massive imputation). Note: For two-level imputation codes
'2' and '-2' are also allowed.
}

\item{where}{A data frame or matrix with logicals of the same dimensions 
as \code{data} indicating where in the data the imputations should be 
created. The default, \code{where = is.na(data)}, specifies that the
missing data should be imputed. The \code{where} argument may be used to 
overimpute observed data, or to skip imputations for selected missing values.}


  \item{visitSequence}{
A vector of integers of arbitrary length, specifying the
column indices of the visiting sequence. The visiting sequence is the column
order that is used to impute the data during one pass through the data. A
column may be visited more than once. All incomplete columns that are used as
predictors should be visited, or else the function will stop with an error.
The default sequence \code{1:ncol(data)} implies that columns are imputed
from left to right. It is possible to specify one of the keywords
\code{'roman'} (left to right), \code{'arabic'} (right to left),
\code{'monotone'} (sorted in increasing amount of missingness) and
\code{'revmonotone'} (reverse of monotone). The keyword should be supplied as
a string and may be abbreviated.
}

\item{blots}{A named \code{list} of \code{alist}'s that can be used 
to pass down arguments to lower level imputation function. The entries
of element \code{blots[[blockname]]} are passed down to the function
called for block \code{blockname}.}


  \item{post}{
A vector of strings with length \code{ncol(data)}, specifying
expressions. Each string is parsed and executed within the \code{sampler()}
function to postprocess imputed values.  The default is to do nothing,
indicated by a vector of empty strings \code{''}.
}

\item{blocks}{List of vectors with variable names per block. List elements 
may be named to identify blocks. Variables within a block are 
imputed by a multivariate imputation method
(see \code{method} argument). By default each variable is placed 
into its own block, which is effectively
fully conditional specification (FCS) by univariate models 
(variable-by-variable imputation). Only variables whose names appear in 
\code{blocks} are imputed. The relevant columns in the \code{where} 
matrix are set to \code{FALSE} of variables that are not block members. 
A variable may appear in multiple blocks. In that case, it is 
effectively re-imputed each time that it is visited.}

\item{formulas}{A named list of formula's, or expressions that
can be converted into formula's by \code{as.formula}. List elements
correspond to blocks. The block to which the list element applies is 
identified by its name, so list names must correspond to block names.
The \code{formulas} argument is an alternative to the 
\code{predictorMatrix} argument that allows for more flexibility in 
specifying imputation models, e.g., for specifying interaction terms.}


  \item{defaultMethod}{
A vector of three strings containing the default
imputation methods for numerical columns, factor columns with 2 levels, and
columns with (unordered or ordered) factors with more than two levels,
respectively. If nothing is specified, the following defaults will be used:
\code{pmm}, predictive mean matching (numeric data) \code{logreg}, logistic
regression imputation (binary data, factor with 2 levels) \code{polyreg},
polytomous regression imputation for unordered categorical data (factor >= 2
levels) \code{polr}, proportional odds model for (ordered, >= 2 levels)
}
  \item{maxit}{
A scalar giving the number of iterations. The default is 5.
}


  \item{seed}{An integer that is used as argument by the \code{set.seed()} for
offsetting the random number generator. Default is to leave the random number
generator alone.}



\item{data.init}{A data frame of the same size and type as \code{data},
without missing data, used to initialize imputations before the start of the
iterative process.  The default \code{NULL} implies that starting imputation
are created by a simple random draw from the data. Note that specification of
\code{data.init} will start the \code{m} Gibbs sampling streams from the same
imputations.}

\item{nnodes}{A scalar indicating the number of nodes for parallel calculation. Default value is 5.}

\item{path.outfile}{A vector of strings indicating the path for redirection of print messages. Default value is NULL, meaning that silent imputation is performed. Otherwise, print messages are saved in the files path.outfile/output.txt. One file per node is generated.}

\item{...}{Named arguments that are passed down to the elementary imputation functions.}
}

\details{
Performs multiple imputation of \code{m} tables in parallel by generating \code{m} seeds, and then by performing multiple imputation by chained equations in parallel from each one. The output is the same as the \code{mice} function of the mice package.
}

\value{
Returns an S3 object of class \code{\link[=mids-class]{mids}} (multiply imputed data set)
}

\references{
Van Buuren, S., Groothuis-Oudshoorn, K. (2011). \code{mice}:
Multivariate Imputation by Chained Equations in \code{R}. \emph{Journal of
Statistical Software}, \bold{45}(3), 1-67.
\url{https://www.jstatsoft.org/article/view/v045i03} <doi:10.18637/jss.v045.i03>

van Buuren, S. (2012). \emph{Flexible Imputation of Missing Data.} Boca
Raton, FL: Chapman & Hall/CRC Press.

Van Buuren, S., Brand, J.P.L., Groothuis-Oudshoorn C.G.M., Rubin, D.B. (2006)
Fully conditional specification in multivariate imputation.  \emph{Journal of
Statistical Computation and Simulation}, \bold{76}, 12, 1049--1064. <doi:10.1080/10629360600810434>

Van Buuren, S. (2007) Multiple imputation of discrete and continuous data by
fully conditional specification.  \emph{Statistical Methods in Medical
Research}, \bold{16}, 3, 219--242. <doi:10.1177/0962280206074463>

Van Buuren, S., Boshuizen, H.C., Knook, D.L. (1999) Multiple imputation of
missing blood pressure covariates in survival analysis.  \emph{Statistics in
Medicine}, \bold{18}, 681--694. <doi:10.1002/(SICI)1097-0258(19990330)18:6<681::AID-SIM71>3.0.CO;2-R>

Brand, J.P.L. (1999) \emph{Development, implementation and evaluation of
multiple imputation strategies for the statistical analysis of incomplete
data sets.} Dissertation. Rotterdam: Erasmus University.
}
\author{
Vincent Audigier \email{vincent.audigier@cnam.fr}
}


\seealso{
\code{\link{mice},\link{parallel}}
}
\examples{

  ##############
  # nhanes (one level data)
  ##############
  data(nhanes, package = "mice")
  #imp <- mice.par(nhanes)
  #fit <- with(data = imp, exp = lm(bmi ~ hyp + chl))
  #summary(pool(fit))
  
  ##############
  #CHEM97Na (Two levels data with 1681 observations and 5 variables)
  ##############
  
  data(CHEM97Na)
  
  ind.clust<-1#index for the cluster variable
  
  #initialisation of the argument predictorMatrix
  predictor.matrix<-mice(CHEM97Na,m=1,maxit=0)$pred
  predictor.matrix[ind.clust,ind.clust]<-0
  predictor.matrix[-ind.clust,ind.clust]<- -2
  predictor.matrix[predictor.matrix==1]<-2
  
  #initialisation of the argument method
  method<-find.defaultMethod(CHEM97Na,ind.clust)
  
  #multiple imputation by chained equations (parallel calculation) [1 minute]
  #(the imputation process can be followed by opening output.txt files in the working directory)
  #res.mice<-mice.par(CHEM97Na,
  #                  predictorMatrix = predictor.matrix,
  #                  method=method,
  #                  path.outfile=getwd())

  
  #multiple imputation by chained equations (without parallel calculation) [4.8 minutes]
  #res.mice<-mice(CHEM97Na,
  #                  predictorMatrix = predictor.matrix,
  #                  method=method)

  
  
  ############
  #IPDNa (Two levels data with 11685 observations and 10 variables)
  ############
  
  data(IPDNa)
  
  ind.clust<-1#index for the cluster variable

  #initialisation of the argument predictorMatrix
  predictor.matrix<-mice(IPDNa,m=1,maxit=0)$pred
  predictor.matrix[ind.clust,ind.clust]<-0
  predictor.matrix[-ind.clust,ind.clust]<- -2
  predictor.matrix[predictor.matrix==1]<-2

  #initialisation of the argument method
  method<-find.defaultMethod(IPDNa,ind.clust)

  #multiple imputation by chained equations (parallel calculation)

  #res.mice<-mice.par(IPDNa,
  #                 predictorMatrix = predictor.matrix,
  #                 method=method,
  #                 path.outfile=getwd())
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{mice}
