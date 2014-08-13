\name{feasible}
\alias{feasible}
\keyword{survey}
\concept{Bayesian survey sampling}
\concept{Polya posterior}
\title{Feasible solution for a probability distribution which must
  satisfy a system of linear equality and inequality constraints.}
\description{This function finds a feasible solution, p=(p1,...,pn), in
  the n-dimensional   simplex of
  probability distributions which must satisfy A1p = b1, A2p <= b2 and
  A3p >= b3. All the components of the bi's must be nonnegative
  In addition  each probability in the solution must
  be at least as big as eps, a small positive number.}
\usage{feasible(A1,A2,A3,b1,b2,b3,eps)}
\arguments{
 \item{A1}{The matrix for the equality constraints.This must always
  contain the constraint that the sum of the pi's is one.}
\item{A2}{The matrix for the <= inequality constraints. This must always
  contain the constraints -pi <= 0, i.e. that the pi's must be
  nonnegative.}
\item{A3}{The matrix for the >= inequality constraints. If there are no
  such constraints A3 must be set equal to NULL.}
\item{b1}{The rhs vector for A1, each component must be nonnegative.}
\item{b2}{The rhs vector for A2, each component must be nonnegative.}
\item{b3}{The rhs vector for A3, each component must be nonnegative.
  If A3 is NULL then b3 must be NULL.}
\item{eps}{A small positive number. Each member of the solution must
  be at least as large as eps. Care must be taken not to choose a value
  of eps which is too large.}
}
\value{The function returns a vector. If the components of the vector
  are positive then the feasible solution is the vector returned,
  otherwise there is no feasible solution.}
\examples{
A1<-rbind(rep(1,7),1:7)
b1<-c(1,4)
A2<-rbind(c(1,1,1,1,0,0,0),c(.2,.4,.6,.8,1,1.2,1.4))
b2<-c(1,2)
A3<-rbind(c(1,3,5,7,9,10,11),c(1,1,1,0,0,0,1))
b3<-c(5,.5)
eps<-1/100
feasible(A1,A2,A3,b1,b2,b3,eps)
}

