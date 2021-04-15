\name{pos.part}
\alias{pos.part}
\title{positive part of a variable}
\description{This is a  function computing the positive part of a variable.}
\usage{
pos.part(x)
}
\arguments{
  \item{x}{A vector of response.}
}
\value{
\item{ps}{The positive part of a variable.}
}
\examples{

pt1 <- proc.time()
##Simulated data
set.seed(12345)
n=200
x = rnorm(n)
obj <- pos.part(x)

proc.time() - pt1
}
