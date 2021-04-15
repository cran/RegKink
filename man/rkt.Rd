\name{rkt}
\alias{rkt}
\title{Estimation for regression kink with a time-varying threshold}
\description{This is a  function estimating regression kink with a time-varying threshold.}
\usage{
rkt(y,x,z,q,r01,r02,r11,r12,stp1,stp2)
}
\arguments{
  \item{y}{A vector of response.}
  \item{x}{A vector of regressor}
  \item{z}{A data matrix of control variables}
  \item{q}{A vector of variable affecting threshold}
  \item{r01}{Lower bounder of parameter space for r0}
  \item{r02}{Upper bounder of parameter space for r0}
  \item{r11}{Lower bounder of parameter space for r1}
  \item{r12}{Upper bounder of parameter space for r1}
  \item{stp1}{Step used in grid search of r0}
  \item{stp2}{Step used in grid search of r1}
}
\value{
A list with the elements
\item{bols}{The OLS estimates when a kink effect is ignored.}
\item{bt}{The regression coefficients when a kink effect is included in the model.}
\item{gammahat0}{The estimated threshold of the constant one in threshold parameters.}
\item{gammahat1}{The estimated threshold of the slop in threshold parameters.}
\item{sig}{The sum of squred errors of the kink model.}
}

\examples{

sta <- proc.time()
##Simulated data
set.seed(12345)
n=200
x = rnorm(n)
q = rnorm(n)
rt = 0.2 - 0.5*q
z = rnorm(n)
x1 = cbind(neg.part(x-rt),pos.part(x-rt),z)
b0 =c(1,2,1)
y = b0[1]*x1[,1]+b0[2]*x1[,2]+b0[3]*x1[,3]+ rnorm(n)



# set grid search paramaters
r01 = 0
r02 = 2
stp1 = 0.1
r11 = -10
r12 = 5
stp2 = 0.1
# estimate the model with a state-dependent threshold
est1 <- rkt(y,x,z,q,r01,r02,r11,r12,stp1,stp2)

proc.time() - sta
}
