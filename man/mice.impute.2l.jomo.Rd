\name{mice.impute.2l.jomo}
\alias{mice.impute.2l.jomo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imputation of univariate missing data by a Bayesian multivariate generalized model based on conjugate priors
}
\description{
Univariate imputation by a Bayesian multivariate generalized model based on conjugate priors. Can be used for a continuous or binary incomplete variable. For continuous variables, the modelling assumes heteroscedasticity for errors. For a binary variable, a probit link and a latent variables framework are used. The method should be used for a variable with large number of clusters and a large number of individuals per cluster.
}
\usage{
mice.impute.2l.jomo(y, ry, x, type, nburn = 200, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Incomplete data vector of length \code{n}
}
  \item{ry}{
Vector of missing data pattern \code{(FALSE=missing, TRUE=observed)}
}
  \item{x}{
Matrix \code{(n x p)} of complete covariates.
}
  \item{type}{
Vector of length \code{ncol(x)} identifying random and class variables. Random variables are identified by a '2'. The class variable (only one is allowed) is coded as '-2'. Random variables also include the fixed effect.
}
  \item{nburn}{
A scalar indicating the number of iterations for the Gibbs sampler. Default is \code{nburn=200}
}
\item{\dots}{
Other named arguments.
}
}
\details{
Contrary to the approach developped in the R jomo package, the imputation is here sequentiallyperformed through a FCS approach, instead of imputing all variables simulatenously. The motivation for such a method is that jomo presents some advantages over other imputation methods, but not always for any type of variables (binary or continuous). By proposing a FCS version of jomo, we allow imputation of mixed variables (continuous and binary), while taking the best of jomo and of other imputation methods. To impute one variable according to this method, other variables are assumed to be full, like in any FCS approach. The imputation function is a direct use of the R function jomo1ran from the jomo package. The argument \code{meth} is tuned to "random" to allow covariance matrices drawn from an inverse Wishart distribution. Only intercept are considered in covariates (\code{X=NULL} and \code{Z=NULL}), while the multivariate outcome corresponds to all variables of the datasets.
}
\value{
A vector of length \code{nmis} with imputations.
}

\author{
Vincent Audigier \email{vincent.audigier@univ-paris-diderot.fr} from the R code of Matteo Quartagno.
}

\references{
Carpenter J.R., Kenward M.G., (2013), Multiple Imputation and its Application. Chapters 3-5-9, Wiley, ISBN: 978-0-470-74052-1.

Yucel R.M., (2011), Random-covariances and mixed-effects models for imputing multivariate multilevel continuous data, Statistical Modelling, 11 (4), 351-370, <doi:10.1177/1471082X100110040>.
}

\seealso{
\code{\link{mice},\link{jomo1ran}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\keyword{mice}
