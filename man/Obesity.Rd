\name{Obesity}
\docType{data}
\alias{Obesity}

\title{A two-level incomplete dataset based on an online obesity survey}
\description{

This synthetic dataset was generated from an online survey on obesity, which collected information on the dietary behavior of 2111 participants. We made the assumption that the data was gathered from five distinct locations or clusters. To account for potential selection bias in the responses related to weight, we simulated the values and observability of this variable using the Heckman selection model within a hierarchical structure.

Additionally, we assumed that in one of the locations, the weight variable was systematically missing. We also introduced missing values for some other variables in the dataset using a Missing at Random (MAR) mechanism.
}


\format{
A dataframe with 2111 observations with the following variables:
 \tabular{rll}{
     \tab Gender \tab a factor variable with two levels: 1 ("Female"), 0 ("Male").\cr
     \tab Age \tab a numeric variable indicating the subject's age in years.\cr
     \tab Height\tab a numeric value with Height in meters.\cr
     \tab FamOb\tab a factor variable describing the subject's family history of obesity with two levels: 1("Yes"), 0("No").\cr
     \tab Weight\tab a numeric variable indicating the subject's weight in kilograms.\cr
     \tab Time\tab a numeric variable indicating the time taken by the subject to respond to the surveys questions in minutes.\cr
     \tab BMI\tab a numeric variable with the subject's body mass index.\cr
     \tab Cluster\tab a numeric variable indexing the cluster.\cr
  }
}


\source{Synthetic data based on the data retrieved from "https://www.kaggle.com/datasets/fabinmndez/obesitydata/"}
\references{
Palechor, F. M., & de la Hoz Manotas, A. (2019). Dataset for estimation of obesity levels based on eating habits and physical condition in individuals from Colombia, Peru and Mexico. Data in brief, 25, 104344.
}

\details{
Data generation code availble on https://github.com/johamunoz/Statsmed_Heckman/blob/main/4.Codes/gendata_Obesity.R
}

\examples{
library(mice)
library(ggplot2)
library(data.table)

data(Obesity)
summary(Obesity)
md.pattern(Obesity)

# Missingness per region (Weight)
dataNA <- setDT(Obesity)[, .(nNA = sum(is.na(Weight)),n=.N), by = Cluster]
dataNA[, propNA:=nNA/n]
dataNA

# Density per region (Weight)
Obesity$Cluster <- as.factor(Obesity$Cluster)
ggplot(Obesity, aes(x = Weight, group=Cluster)) +
  geom_histogram(aes(color = Cluster,fill= Cluster),
                 position = "identity", bins = 30) + 
                 facet_grid(Cluster~.)
}

\keyword{datasets}