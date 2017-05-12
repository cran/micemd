\name{mice.impute.2l.2stage.norm}
\alias{mice.impute.2l.2stage.norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imputation by a two-level heteroscedastic normal model based on a two-stage estimator
}
\description{
Imputes univariate two-level continuous variable from a heteroscedastic normal model. The imputation method is based on a two-stage estimator: at step 1, a linear regression model is fitted to each observed cluster; at step 2, estimates obtained from each cluster are combined according to a linear random effect model.

}
\usage{
mice.impute.2l.2stage.norm(y, ry, x, type, method_est = "mm", ...)
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
  \item{\dots}{
Other named arguments.
}
}
\details{
Imputes univariate two-level continuous variable from a heteroscedastic normal model. The imputation method is based on a two-stage estimator: at step 1, a linear regression model is fitted to each observed cluster; at step 2, estimates obtained from each cluster are combined according to a linear random effect model. Two possibilities are available to combine estimates at stage 2: by default, parameters of the linear random effect model are estimated according to the method of moments (MM), otherwise, parameters of the linear random effect model can be estimated according to the restricted maximum likelihood estimator (REML). The variability on the parameters of the imputation is propagated according to an asymptotic strategy requiring notably a large number of clusters. Compared to the REML version, the MM version is quicker to perform, but it provides less theoretical garanties. Nevertheless, simulations studies show that both versions lead to similar inferences (Resche-Rigon, M. and White, I. R. (2016)).
}
\value{
A vector of length \code{nmis} with imputations.
}

\references{
Resche-Rigon, M. and White, I. R. (2016). Multiple imputation by chained equations for systematically
and sporadically missing multilevel data. Statistical Methods in Medical Research.
To appear. <doi:10.1177/0962280216666564>

Audigier, V., White, I. ,  Jolani ,S.  Debray, T.,  Quartagno, M., Carpenter, J., van Buuren, S. and Resche-Rigon, M.
Multiple imputation for multilevel data with continuous and binary variables. ArXiv e-prints <arXiv:1702.00971>, 2017.

}

\author{
Vincent Audigier \email{vincent.audigier@univ-paris-diderot.fr}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{
\code{\link{mice},\link{mice.impute.2l.2stage.pmm},\link{mice.impute.2l.glm.norm},\link{mice.impute.2l.jomo}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{mice}
