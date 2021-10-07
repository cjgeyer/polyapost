\name{feasible}
\alias{feasible}
\keyword{survey}
\concept{Bayesian survey sampling}
\concept{Polya posterior}
\title{Feasible Solution for a Probability Distribution which must
  Satisfy a System of Linear Equality and Inequality Constraints.}
\description{This function finds a feasible solution,
  \eqn{p=(p_1,\ldots,p_n)}{p=(p1,...,pn)}, in the \eqn{n}-dimensional simplex of
  probability distributions which must satisfy
  \eqn{A_1 p = b_1}{A1 p = b1},
  \eqn{A_2 p = b_2}{A2 p = b2}, and
  \eqn{A_3 p = b_3}{A3 p = b3},
  All the components of the \eqn{b_i}{bi} must be nonnegative
  In addition  each probability in the solution must
  be at least as big as \code{eps}, a small positive number.}
\usage{feasible(A1,A2,A3,b1,b2,b3,eps)}
\arguments{
 \item{A1}{The matrix for the equality constraints.This must always
  contain the constraint \code{sum(p) == 1}.}
\item{A2}{The matrix for the \code{<=} inequality constraints.  This must always
  contain the constraints \code{-p <= 0}.}
\item{A3}{The matrix for the \code{>=} inequality constraints.  If there are no
  such constraints \code{A3} must be set equal to \code{NULL}.}
\item{b1}{The rhs vector for \code{A1}, each component must be nonnegative.}
\item{b2}{The rhs vector for \code{A2}, each component must be nonnegative.}
\item{b3}{The rhs vector for \code{A3}, each component must be nonnegative.
  If \code{A3} is \code{NULL} then \code{b3} must be \code{NULL}.}
\item{eps}{A small positive number.  Each member of the solution must
  be at least as large as \code{eps}.  Care must be taken not to choose a value
  of \code{eps} which is too large.}
}
\value{The function returns a vector.  If the components of the vector
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

