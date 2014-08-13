\name{wtpolyap}
\alias{wtpolyap}
\keyword{survey}
\concept{Bayesian survey sampling}
\concept{Polya posterior}
\title{Polya sampling from an urn with possibly unequal weights}
\description{Consider an urn containing
  a finite set of values along with their respective positive weights.
  An item is selected at random from the urn with probability
  proportional to its weight. Then it is returned to the urn and its
  weight is increased by one. The process is repeated on the
  adjusted urn. We continue until the total weight in the urn has
  been increased by k. 
  The original composition of the urn along with the k selected
  values, in order, are returned.}
\usage{wtpolyap(ysamp, wts, k)}
\arguments{
  \item{ysamp}{A vector of real numbers which make up the urn.}
  \item{wts}{A vector of positive weights which defines the initial
    probability of selection.}
  \item{k}{A positive integer which specifies the number of Polya
   samples taken from the urn where after each draw the weight of the selected
   item is increased by one.}
}
\value{The returned value is a vector of length equal to the length of
  the sample plus k}
\examples{
wtpolyap(c(0,1,2),c(0.5,1,1.5),22)
}
