\name{IPDNa}
\alias{IPDNa}
\docType{data}
\title{
A simulated Individual Patient Data (IPD) meta-analysis with missing values.
}
\description{
This dataset is a simulated version of an IPD meta-analysis consisting of 28 studies focusing on risk factors in acute heart failure (GREAT, 2013). Each study includes a list of patient characteristics and potential risk factors. Each of them is incomplete, leading to sporadically missing values (Resche-Rigon, et al 2013). In addition, some variables have been collected on some studies only, leading to systematically missing values. More details on the original dataset are provided in Audigier et al. (2018). To mimic the real data, a general location model has been fitted on each study (Schafer, 1997). Then, each study has been generated according to the estimated parameters. Finally, missing values have been allocated similarly to the original dataset.
}
\usage{data("IPDNa")}
\format{
  A data frame with 11685 observations on the following 10 variables.
  \describe{
    \item{\code{centre}}{a numeric indexing the center where the study is conducted}
    \item{\code{gender}}{a factor with levels \code{0} \code{1}}
    \item{\code{bmi}}{a numeric vector indicating the body mass index}
    \item{\code{age}}{a numeric vector indicating the age}
    \item{\code{sbp}}{a numeric vector indicating the systolic blood pressure}
    \item{\code{dbp}}{a numeric vector indicating the diastolic blood pressure}
    \item{\code{hr}}{a numeric vector indicating the heart rate}
    \item{\code{lvef}}{a numeric vector indicating the ventricular ejection fraction}
    \item{\code{bnp}}{a numeric vector indicating the level of the brain natriuretic peptide biomarker}
    \item{\code{afib}}{a factor with levels \code{0} \code{1} indicating the atrial fibrillation}
  }
}
\details{
For more details, see Audigier et al. (2018)
}
\source{
GREAT Network (2013).  Managing acute heart failure in the ed - case studies
from the acute heart failure academy.  http://www.greatnetwork.org
}
\references{
Audigier, V., White, I. ,  Jolani ,S.  Debray, T.,  Quartagno, M., Carpenter, J., van Buuren, S. and Resche-Rigon, M. Multiple imputation for multilevel data with continuous and binary variables (2018). Statistical Science. \doi{10.1214/18-STS646}.

Resche-Rigon, M., White, I. R., Bartlett, J., Peters, S., Thompson, S., and on behalf of the PROG-IMT Study Group (2013). Multiple imputation for handling systematically missing confounders in meta-analysis of individual participant data. Statistics in Medicine, 32(28):4890-4905. \doi{10.1002/sim.5894}

Schafer, J. L. (1997) Analysis of Incomplete Multivariate Data. Chapman & Hall, Chapter 9.
}
\examples{
data(IPDNa)

#summary
summary(IPDNa)

#summary per study
by(IPDNa, IPDNa$centre, summary)
}
\keyword{datasets}
