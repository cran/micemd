\name{mice.impute.2l.2stage.pmm}
\alias{mice.impute.2l.2stage.pmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictive mean matching imputation for two-level variable
}
\description{
Similarly to mice.impute.2l.stage.norm, this function imputes univariate two-level continuous variable from a heteroscedastic normal model. The difference consists in replacing missing values by observed values instead of adding a parametric noise to the prediction of a linear model with random effects (as done in mice.impute.2l.stage.norm.mm and mice.impute.2l.stage.norm.reml).
}
\usage{
mice.impute.2l.2stage.pmm(y, ry, x, type,
                              method_est = "mm",
                              incluster = FALSE,
                              kpmm = 5, ...)
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
  \item{method_est}{
Vector of string given the version of the estimator to used. Choose \code{method_est="reml"} for restricted maximum likelihood estimator or \code{method_est="mm"} for the method of moments. By default \code{method_est="mm"}.
}
  \item{incluster}{
Boolean indicating if the imputed values are drawn from the cluster or from the full dataset. By default imputed values are drawn from all available clusters \code{incluster=FALSE}.
}
  \item{kpmm}{
The size of the donor pool among which a draw is made. The default is \code{k = 5}.
}
\item{\dots}{
Other named arguments.
}
}
\details{
Imputes univariate two-level continuous variable from observed values. The imputation method is based on a two-stage estimator: at step 1, a linear regression model is fitted to each observed cluster; at step 2, estimates obtained from each cluster are combined according to a linear random effect model. To combine estimates at stage 2, parameters of the linear random effect model are estimated according to the method of moments or according to the restricted maximum likelihood estimator. The variability on the parameters of the imputation is propagated according to an asymptotic strategy requiring a large number of clusters. The sample variability is reflected by using a predictive mean matching approach, meaning that missing values are imputed by a draw from observed values. The pool of \code{k} donors is defined according to the Manhattan distance between the prediction of the observation which is imputed and the predictions of other available observations (matching of type 2). The pool can be restricted to the cluster of the individual that is imputed or from all clusters. By drawing values inside the cluster, the heteroscedasticity assumption is preserved. Otherwise, the sample variability of imputed values is the same for all clusters, which strengthen the homoscedasticity assumption. Among the pool of \code{k} donors, the selected one is drawn at random.
}
\value{
Numeric vector of length  \code{sum(!ry)} with imputations
}
\references{
Rubin, D.B. (1987). Multiple imputation for nonresponse in surveys. New York: Wiley.

Resche-Rigon, M. and White, I. R. (2016). Multiple imputation by chained equations for systematically
and sporadically missing multilevel data. Statistical Methods in Medical Research.
To appear. <doi:10.1177/0962280216666564>

Audigier, V., White, I. ,  Jolani ,S.  Debray, T.,  Quartagno, M., Carpenter, J., van Buuren, S. and Resche-Rigon, M.
Multiple imputation for multilevel data with continuous and binary variables (2018). Statistical Science. <doi:10.1214/18-STS646>.

}
\author{
Vincent Audigier \email{vincent.audigier@cnam.fr}
}
\note{
This method is experimental.
}

\seealso{
\code{\link{mice.impute.2l.2stage.norm}}
}
