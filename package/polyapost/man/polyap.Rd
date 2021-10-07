\name{polyap}
\alias{polyap}
\keyword{survey}
\concept{Bayesian survey sampling}
\concept{Polya posterior}
\title{Polya Sampling from an Urn}
\description{Consider an urn containing a finite set of values.
  An item is selected at random from the urn.  Then it is returned
  to the urn along with another item with the same value.  Next
  a value is selected at random from the reconstituted urn
  and it and a copy our returned to the urn.  This process is
  repeated until \eqn{k} additional items have been added to the
  original urn.  The original composition of the urn along with the selected
  values, in order, are returned.}
\usage{polyap(ysamp, k)}
\arguments{
\item{ysamp}{A vector of real numbers which make up the urn.}
\item{k}{A positive integer which specifies the number of items
  added to the original composition of the urn.}
}
\value{The returned value is a vector of length equal to the length of
  \code{ysamp} plus \code{k}.}
\examples{
polyap(c(0,1),20)
}
