\name{mice.impute.2l.2stage.heckman}
\alias{mice.impute.2l.2stage.heckman}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Imputation based on Heckman model for multilevel data.
}
\description{
Imputes outcome and predictor variables that follow an MNAR mechanism according to Heckman's model and come from a multilevel database such as individual participant data with systematically and sporadically missing values.
}
\usage{
mice.impute.2l.2stage.heckman(y, ry, x, wy = NULL, type,
pmm = FALSE, ypmm = NULL, meta_method = "reml", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Vector to be imputed
}
  \item{ry}{
Logical vector of length \code{length(y)} indicating the subset \code{y[ry]} of elements in \code{y} to which the imputation model is fitted. The \code{ry} generally distinguishes the observed (\code{TRUE}) and missing values (\code{FALSE}) in \code{y}.
}
  \item{x}{
Numeric design matrix with \code{length(y)} rows with predictors for \code{y}. Matrix \code{x} may have no missing values.
}
  \item{wy}{
Logical vector of length \code{length(y)}. A \code{TRUE} value indicates locations in \code{y} for which imputations are created.
}
  \item{type}{
Type of the variable in the prediction model {0: No predictor,
 1: Predictor in both the outcome and selection,-2: Cluster id (study id), -3: Predictor only in the selection model, -4: Predictor only in the outcome model.}}
  \item{pmm}{
Predictive mean matching can be applied only for missing continuous variables: "TRUE","FALSE".
}
  \item{ypmm}{
Continuous vector of donor values of y to perform the predictive mean matching, in case ypmm is not provided, the observable values of y are used.
}
  \item{meta_method}{
Meta_analysis estimation method for random effects: "ml" (maximum likelihood), "reml" (restricted maximum likelihood) or "mm" method of moments.
}
  \item{\dots}{
Other named arguments. Not used.
}}
\details{
Imputes systematically and sporadically missing binary and continuous univariate variables that follow a MNAR mechanism according to the Heckman selection model and come from a clustered dataset. The imputation method uses a two-stage approach in which the Heckman model parameters at the cluster level are estimated using the copula method.
}
\value{
Vector with imputed data, of type binary or continuous
}
\references{
Munoz J,Hufstedler H,Gustafson P, Barnighausen T, De Jong V, Debray T. Dealing with missing data using the Heckman selection model: methods primer for epidemiologists.IJE,December 2022. DOI:	10.1093/ije/dyac237.

Munoz J, Egger M, Efthimiou O, Audigier V, De Jong V, Debray T. Multiple imputation of incomplete multilevel data using Heckman selection models, Jan 2023, https://doi.org/10.48550/arXiv.2301.05043.
}
\author{
Julius Center for Health Sciences and Primary Care, University Medical Center Utrecht, Utrecht University, 2022 \email{j.munozavila@umcutrecht.nl}
}
\note{
Missing binary variables should be included as two-level factor type variables in the incomplete dataset. Cluster variable should be included as numeric variable in the dataset, and defined as -2 in the predictor matrix. When the cluster variable is not specified, the imputation method is based on a simple Heckman model, i.e. without taking into account the hierarchical structure. In case the Heckman model cannot be estimated at the cluster level, the imputation method will be based on the simple Heckman model.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mice}}
}
\examples{
require(mice)
data(data_heckman)

#####################
# Initialisation
#####################

# define imputation methods for each incomplete variables
ini <- mice(data_heckman, maxit = 0)
meth <- ini$method 
meth[c("Age","Height")] <- c("2l.norm", "2l.norm") 
meth[c("FAVC")] <- "2l.glm.bin" 
meth["Weight"] <- "2l.2stage.heckman" 

# set type of predictor variable, here Weight variable is assumed an MNAR variable
# as Age and Gender are predictors in both selection and outcome model are set as 1

pred <- ini$predictorMatrix
pred[,"Cluster"] <- -2 # Cluster variable
pred["Weight","Time"]  <- -3 # Variable only affects the selection model (Exclusion restriction)
pred["Weight",c("Height","FAVC")]  <- -4 # Variables only affect the outcome model

#####################
# multiple imputation
#####################
#imp <- mice(data = data_heckman, meth = meth, pred = pred, seed = 123)

#-----------------
# Summary weight
#-----------------

#summary(complete(imp,"long")$Weight)

#-----------------
# Model weight
#-----------------

#library(broom.mixed)
#model_MNAR <- with(imp, lmer(Weight~Gender+Age+Height+FAVC+(1|Cluster)))
#summary(pool(model_MNAR))

#####################
# multiple imputation with pmm
#####################

#imp_pmm <- mice(data = data_heckman, meth = meth, pred = pred,
#pmm=TRUE, ypmm=seq(35,180,0.1), seed = 123)
#model_MNAR_pmm <- with(imp_pmm, lmer(Weight~Gender+Age+Height+FAVC+(1|Cluster)))
#summary(pool(model_MNAR_pmm))
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
