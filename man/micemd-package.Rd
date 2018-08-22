\name{micemd-package}
\alias{micemd-package}
\alias{micemd}
\docType{package}
\title{
Multiple Imputation by Chained Equations with Multilevel Data

}
\description{
Addons for the mice package to perform multiple imputation using chained equations with two-level data. Includes imputation methods specifically handling sporadically and systematically missing values (Resche-Rigon et al. 2013). Imputation of continuous, binary or count variables are available. Following the recommendations of Audigier, V. et al (2018), the choice of the imputation method for each variable can be facilitated by a default choice tuned according to the structure of the incomplete dataset. Allows parallel calculation for mice.
}

\author{
Vincent Audigier, Matthieu Resche-Rigon

Maintainer: Vincent Audigier <vincent.audigier@cnam.fr>
}
\references{
Audigier, V., White, I. ,  Jolani ,S.  Debray, T.,  Quartagno, M., Carpenter, J., van Buuren, S. and Resche-Rigon, M.
Multiple imputation for multilevel data with continuous and binary variables (2018). Statistical Science. <doi:10.1214/18-STS646>.


Jolani, S., Debray, T. P. A., Koffijberg, H., van Buuren, S., and Moons, K. G. M. (2015). Imputation
of systematically missing predictors in an individual participant data meta-analysis: a
generalized approach using MICE. Statistics in Medicine, 34(11):1841{1863}. <doi:10.1002/sim.6451>

Quartagno, M. and Carpenter, J. R. (2016). jomo: A package for Multilevel Joint Modelling
Multiple Imputation.

Quartagno, M. and Carpenter, J. R. (2016). Multiple imputation for IPD meta-analysis: allowing
for heterogeneity and studies with missing covariates. Statistics in Medicine, 35(17):2938{
2954}. <doi:10.1002/sim.6837>

Resche-Rigon, M. and White, I. R. (2016). Multiple imputation by chained equations for systematically
and sporadically missing multilevel data. Statistical Methods in Medical Research.
To appear. <doi:10.1177/0962280216666564>

Resche-Rigon, M., White, I. R., Bartlett, J., Peters, S., Thompson, S., and on behalf of the
PROG-IMT Study Group (2013). Multiple imputation for handling systematically missing
confounders in meta-analysis of individual participant data. Statistics in Medicine,
32(28):4890{4905}. <doi:10.1002/sim.5894>

Van Buuren, S., Groothuis-Oudshoorn, K. (2011). mice: Multivariate Imputation by Chained Equations in R. Journal of Statistical Software, 45(3), 1-67. <doi:10.18637/jss.v045.i03> http://www.jstatsoft.org/v45/i03/
}

\keyword{package}
\keyword{multilevel}
\keyword{imputation}
\keyword{meta-analyse}
\seealso{
\code{\link{mice}}
}
\examples{

  require(lme4)
  data(CHEM97Na)
  
  ind.clust<-1#index for the cluster variable
  
  #initialisation of the argument predictorMatrix
  temp<-mice(CHEM97Na,m=1,maxit=0)
  temp$pred[ind.clust,ind.clust]<-0;temp$pred[-ind.clust,ind.clust]<- -2
  temp$pred[temp$pred==1]<-2
  predictor.matrix<-temp$pred
  
  #initialisation of the argument method
  method<-find.defaultMethod(CHEM97Na,ind.clust)
  
  #multiple imputation by chained equations (parallel calculation) [time consumming]
  #res.mice<-mice.par(CHEM97Na,predictorMatrix = predictor.matrix,method=method)
  
  #check convergence
  #plot(res.mice)
  
  #analysis (apply a generalized linear mixed effects model to each imputed dataset)
  #ana<-with(res.mice,expr=glmer(Score~Sex+GSCE+(1|School),
  #                            family="poisson",
  #                            control=glmerControl(optimizer = "bobyqa")))
                              
  #check the number of generated tables
  #plot(ana)

  #pooling
  #res.pool<-pool(ana)
  #summary(res.pool)
  
}
