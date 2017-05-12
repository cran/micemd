\name{CHEM97Na}
\alias{CHEM97Na}
\docType{data}
\title{
An incomplete two-level dataset which consists of A/AS-level examination data from England
}
\description{
This dataset is an extract of the CHEM97 dataset (Fielding, A. et al, 2003) dealing with point scores of 31,022 pupils grouped in 2,280 schools. CHEM97Na reports point score for Schools with more than 70 pupils only, i.e. 1681 pupils grouped in 18 schools. Systematically missing values and sporadically missing values have been added according to a missing completely at random (MCAR) mechanism (Little R.J.A. and Rubin D.B., 2002). Systematically missing values are values that are missing for all pupils of a same school, while sporadically missing values are values which are missing for an individual only (Resche-Rigon, et al 2013).
}
\usage{data("CHEM97Na")}
\format{
  A data frame with 1681 observations on the following 5 variables.
  \describe{
    \item{\code{School}}{a numeric indexing the School}
    \item{\code{Sex}}{a factor with levels \code{M} \code{F}}
    \item{\code{Age}}{a numeric indicating the age in months}
    \item{\code{GSCE}}{a numeric vector indicating the point score at the General Certificate of Secondary Education}
    \item{\code{Score}}{a numeric vector indicating the point score on A-level Chemistry in 1997}
  }
}
\details{
For more details, see Fielding, A. et al (2003).
}
\source{
Fielding, A., Yang, M., and Goldstein, H.(2003). Multilevel ordinal models for examination grades. Statistical Modelling, 3 (2): 127-153.

Available at http://www.bristol.ac.uk/cmm/learning/mmsoftware/data-rev.html#chem97
}

\references{
Fielding, A., Yang, M., and Goldstein, H. (2003). Multilevel ordinal models for examination grades. Statistical Modelling, 3 (2): 127-153. <doi:10.1191/1471082X03st052oa>

Resche-Rigon, M., White, I. R., Bartlett, J., Peters, S., Thompson, S., and on behalf of the
PROG-IMT Study Group (2013). Multiple imputation for handling systematically missing
confounders in meta-analysis of individual participant data. Statistics in Medicine,
32(28):4890{4905}. <doi:10.1002/sim.5894>

Little R.J.A., Rubin D.B. (2002) Statistical Analysis with Missing Data. Wiley series in probability and statistics, New-York
}
\examples{
data(CHEM97Na)

#summary
summary(CHEM97Na)

#summary per School
by(CHEM97Na,CHEM97Na$School,summary)

}
\seealso{
\code{\link[VIM]{matrixplot}}
}
\keyword{datasets}
