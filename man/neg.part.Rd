\name{neg.part}
\alias{neg.part}
\title{negtive part of a variable}
\description{This is a  function computing the negtive part of a variable.}
\usage{
neg.part(x)
}
\arguments{
  \item{x}{A vector of response.}
}

\value{
\item{ne}{The negative part of a variable.}
}

\examples{

pt1 <- proc.time()
##Simulated data
set.seed(12345)
n=200
x = rnorm(n)
obj <- neg.part(x)

proc.time() - pt1
}
