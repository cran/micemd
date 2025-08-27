\name{mice.impute.2l.glm.bin}
\alias{mice.impute.2l.glm.bin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imputation of univariate missing data using a Bayesian logistic mixed model based on non-informative prior distributions
}
\description{
Imputes univariate missing data using a Bayesian logistic mixed model based on non-informative prior distributions. The method is dedicated to a binary outcome stratified in severals clusters. Should be used with few clusters and few individuals per cluster. Can be very slow to perform otherwise.
}
\usage{
mice.impute.2l.glm.bin(y, ry, x, type, ...)
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
\item{\dots}{
Other named arguments.
}
}

\details{
Imputes univariate missing data using a Bayesian logistic mixed model based on non-informative prior distributions. The variability on the parameters of the imputation is propagated according to an explicit Bayesian modelling. More precisely, improper prior distributions are used for regression coefficients and covariance matrix of random effects. The method is recommended for datasets with a small number of clusters and a small number of individuals per cluster. Otherwise, the method can be very slow to perform.
}
\value{
A vector of length \code{nmis} with imputations.
}

\references{
Jolani, S., Debray, T. P. A., Koffijberg, H., van Buuren, S., and Moons, K. G. M. (2015). Imputation of systematically missing predictors in an individual participant data meta-analysis: a generalized approach using MICE. Statistics in Medicine, 34(11):1841-1863. \doi{10.1002/sim.6451}

Audigier, V., White, I. ,  Jolani ,S.  Debray, T.,  Quartagno, M., Carpenter, J., van Buuren, S. and Resche-Rigon, M. Multiple imputation for multilevel data with continuous and binary variables (2018). Statistical Science. \doi{10.1214/18-STS646}.
}
\author{
Vincent Audigier \email{vincent.audigier@cnam.fr} from the R code of Shahab Jolani.
}
\seealso{
\code{\link[mice]{mice}} \code{\link{mice.impute.2l.2stage.bin}} \code{\link{mice.impute.2l.jomo}}
}

\keyword{mice}
