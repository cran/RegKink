\name{reg}
\alias{reg}
\title{linear regression}
\description{This is a  function computing OLS estimates of linear model.}
\usage{
reg(X,y)
}
\arguments{
  \item{y}{A vector of response.}
  \item{X}{data matrix}
}
\value{
\item{bols}{The ols estimates of a linear regression model.}
}
\examples{

pt1 <- proc.time()
##Simulated data
set.seed(12345)
n=200
e = rnorm(n)
X = cbind(1,rnorm(n))
b0 =c(1,2)
y = X%*%b0 + e
obj <- reg(X,y)

proc.time() - pt1
}
