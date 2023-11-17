\name{mice.impute.2l.2stage.heckman}
\alias{mice.impute.2l.2stage.heckman}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imputation based on Heckman model for multilevel data.
}
\description{
Imputes both outcome or predictor incomplete variables that follow an indirectly non-ignorable Missing Not at Random (MNAR) mechanism, i.e., the likelihood of a missing value in the incomplete variable depends on other unobserved variable(s) that are also correlated with the incomplete variable. This imputation is based on Heckman's selection model and is suitable for multilevel databases, such as individual participant data, with both systematic or sporadic missing data.}
\usage{
mice.impute.2l.2stage.heckman(y, ry, x, wy = NULL, type,
pmm = FALSE, ypmm = NULL, meta_method = "reml", pred_std = FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Vector to be imputed
}
  \item{ry}{
A logical vector of length \code{length(y)} indicating the subset \code{y[ry]} of elements in \code{y} to which the imputation model is fitted. The \code{ry} generally distinguishes the observed (\code{TRUE}) and missing values (\code{FALSE}) in \code{y}.
}
  \item{x}{
A numeric design matrix with \code{length(y)} rows with predictors for \code{y}. Matrix \code{x} may have no missing values.
}
  \item{wy}{
A logical vector of length \code{length(y)}. A \code{TRUE} value indicates locations in \code{y} for which imputations are created.
}
  \item{type}{
Type of the variable in the prediction model, which can be one of the following: No predictor (0), Cluster variable (-2), Predictor in both the outcome and selection equation (2), Predictor only in the selection equation (-3), Predictor only in the outcome equation (-4). In this method all predictors are considered random variables that also included the fixed effect.}
  \item{pmm}{
A logical value that specifies whether the predictive mean matching method is applied.(default = "FALSE"). This method is only applicable to missing continuous variables.
}
  \item{ypmm}{
A continuous vector of donor values for y used in the predictive mean matching method. if ypmm is not provided, the observable values of y are used as donors.
}
  \item{meta_method}{
A character value that indicates the method for estimating meta_analysis random effects: "ml" (maximum likelihood), "reml" (restricted maximum likelihood) or "mm" (method of moments).
}
  \item{pred_std}{
A logical value that indicates whether internally standardize the set of predictor variables (default = FALSE).
}
  \item{\dots}{
Other named arguments. Not used.
}}
\details{
This function imputes systematically and sporadically missing binary and continuous univariate variables that follow an MNAR mechanism according to the Heckman selection model. It is specifically designed for clustered datasets. The imputation method employs a two-stage approach in which the Heckman model parameters at the cluster level are estimated using the copula method.
}
\value{
Vector with imputed data, of type binary or continuous type
}
\references{
Munoz J,Hufstedler H,Gustafson P, Barnighausen T, De Jong V, Debray T. Dealing with missing data using the Heckman selection model: methods primer for epidemiologists.IJE,December 2022. \doi{10.1093/ije/dyac237}.

Munoz J, Egger M, Efthimiou O, Audigier V, De Jong V, Debray T. Multiple imputation of incomplete multilevel data using Heckman selection models, Jan 2023, \doi{10.48550/arXiv.2301.05043}.
}
\author{
Julius Center Methods Group UMC, 2022 \email{J.MunozAvila@umcutrecht.nl}
}
\note{
Missing binary variables should be included as two-level factor type variables in the incomplete dataset.The cluster variable should be included as a numeric variable in the dataset. When the cluster variable is not specified, the imputation method defaults to a simple Heckman model, which does not take in account the hierarchical structure. In cases where the Heckman model cannot be estimated at the hierarchical level, the imputation method reverts to the simple Heckman model.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mice}}
}
\examples{
require(mice)
require(nlme)
require(broom.mixed)
require(parallel)

# Load dataset
data(Obesity) 

# Define imputation methods for each incomplete variables
meth <- find.defaultMethod(Obesity, ind.clust = 1)

# Modify some of the proposed imputation methods
# Deterministic imputation
meth["BMI"] <- "~ I(Weight / (Height)^2)" 
meth["Age"] <- "2l.2stage.pmm" 

# Set method, here Weight variable is assumed an MNAR variable
# Weight imputed with the Heckman method
meth["Weight"] <- "2l.2stage.heckman" 

# Set type of predictor variable, 
# All covariates are included in both outcome and selection equation
ini <- mice(Obesity, maxit = 0)
pred <- ini$pred
pred[,"Time"] <- 0
pred[,"Cluster"] <- -2
pred[pred == 1] <- 2

# Time was used as exclusion restriction variable 
pred["Weight","Time"]  <- -3  
# Deterministic imputation, to avoid circular predictions
pred[c("Height", "Weight"), "BMI"] <- 0 

# Imputation of continuous variables (time consumming)

# nnodes <- detectCores()
# imp <- mice.par(Obesity, meth = meth, pred = pred, m=10, seed = 123,
#                 nnodes = nnodes)
# summary(complete(imp,"long")$Weight)

# Imputation of continuous variables using the predictor mean matching method.
# Imputed values fall within the range of observable variables.

# imp_pmm <- mice.par(Obesity, meth = meth, pred = pred, m = 10,
#                     seed = 123, pmm=TRUE, nnodes = nnodes) 
# summary(complete(imp_pmm,"long")$Weight)

# Fit the model

# model_MNAR <- with(imp,lme( BMI ~ Age + FamOb + Gender,random=~1+Age|Cluster))
# model_MNAR_pmm <- with(imp_pmm,lme( BMI ~ Age + FamOb + Gender,random=~1+Age|Cluster))

# summary(pool(model_MNAR))
# summary(pool(model_MNAR_pmm))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
