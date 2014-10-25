\name{hitrun}
\alias{hitrun}
\alias{hitrun.default}
\alias{hitrun.hitrun}
\title{Hit and Run Algorithm for Constrained Dirichlet Distribution}
\description{
    Markov chain Monte Carlo for equality and inequality constrained
    Dirichlet distribution using a hit and run algorithm.
}
\usage{
hitrun(alpha, \dots)

\S3method{hitrun}{default}(alpha, a1 = NULL, b1 = NULL, a2 = NULL, b2 = NULL,
     nbatch = 1, blen = 1, nspac = 1, outmat = NULL, debug = FALSE, \dots)

\S3method{hitrun}{hitrun}(alpha, nbatch, blen, nspac, outmat, debug, \dots)
}
\arguments{
  \item{alpha}{parameter vector for Dirichlet distribution.  Alternatively,
      an object of class \code{"hitrun"} that is the result of a previous
      invocation of this function, in which case this run continues where
      the other left off.}
  \item{nbatch}{the number of batches.}
  \item{blen}{the length of batches.}
  \item{nspac}{the spacing of iterations that contribute to batches.}
  \item{a1}{a numeric or character matrix or \code{NULL}. See details.}
  \item{b1}{a numeric or character vector or \code{NULL}. See details.}
  \item{a2}{a numeric or character matrix or \code{NULL}. See details.}
  \item{b2}{a numeric or character vector or \code{NULL}. See details.}
  \item{outmat}{a numeric matrix, which controls the output.  If \code{p}
      is the constrained Dirichlet random vector being simulated, then
      \code{outmat \%*\% p} is the functional of the state that is averaged.
      May be \code{NULL}, in which case the identity matrix is used.}
  \item{debug}{if \code{TRUE}, then additional output useful for debugging
      is produced.}
  \item{\dots}{ignored arguments.  Allows the two methods to have different
      arguments.  You cannot change the Dirichlet parameter or the constraints
      (hence cannot change the target distribution) when using the method
      for class \code{"hitrun"}.}
}
\details{
Runs a hit and run algorithm (for which see the references)
producing a Markov chain with equilibrium distribution having a Dirichlet
distribution with parameter vector \code{alpha} constrained to lie in the
subset of the unit simplex consisting of \code{x} satisfying
\preformatted{
    a1 \%*\% x <= b1
    a2 \%*\% x == b2
}
Hence if \code{a1} is \code{NULL} then so must be \code{b1}, and vice versa,
and similarly for \code{a2} and \code{b2}.

If any of \code{a1}, \code{b1}, \code{a2}, \code{b2} are of type
\code{"character"}, then they must be valid GMP (GNU multiple precision)
rational, that is, if run through \code{\link{q2q}}, they do not
give an error.  This allows constraints to be represented exactly
(using infinite precision rational arithmetic) if so desired.
}
\value{
  an object of class \code{"hitrun"},
  which is a list containing at least the following components:
  \item{batch}{\code{nbatch} by \code{p} matrix, the batch means, where
      \code{p} is the row dimension of \code{outmat}.}
  \item{initial}{initial state of Markov chain.}
  \item{final}{final state of Markov chain.}
  \item{initial.seed}{value of \code{.Random.seed} before the run.}
  \item{final.seed}{value of \code{.Random.seed} after the run.}
  \item{time}{running time from \code{system.time()}.}
  \item{alpha}{the Dirichlet parameter vector.}
  \item{nbatch}{the argument \code{nbatch} or \code{obj$nbatch}.}
  \item{blen}{the argument \code{blen} or \code{obj$blen}.}
  \item{nspac}{the argument \code{nspac} or \code{obj$nspac}.}
  \item{outmat}{the argument \code{outmat} or \code{obj$outmat}.}
%   Description of additional output when \code{debug = TRUE} can be
%   found in the vignette \code{debug} (\url{../doc/debug.pdf}).
}
\references{
Belisle, C. J. P., Romeijn, H. E. and Smith, R. L. (1993)
Hit-and-run algorithms for generating multivariate distributions.
\emph{Mathematics of Operations Research}, \bold{18}, 255--266.

Chen, M. H. and Schmeiser, B. (1993)
Performance of the Gibbs, hit-and-run, and Metropolis samplers.
\emph{Journal of Computational and Graphical Statistics}, \bold{2}, 251--272.

}
\seealso{\code{\link{ConvertGMP}} and \code{\link{ArithmeticGMP}}}
\examples{
# Bayesian inference for discrete probability distribution on {1, ..., d}
# state is probability vector p of length d
d <- 10
x <- 1:d
# equality constraints
#     mean equal to (d + 1) / 2, that is, sum(x * p) = (d + 1) / 2
# inequality constraints
#     median less than or equal to (d + 1) / 2, that is,
#         sum(p[x <= (d + 1) / 2]) <= 1 / 2
a2 <- rbind(x)
b2 <- (d + 1) / 2
a1 <- rbind(as.numeric(x <= (d + 1) / 2))
b1 <- 1 / 2
# simulate prior, which Dirichlet(alpha)
# posterior would be another Dirichlet with n + alpha - 1,
#    where n is count of IID data for each value
alpha <- rep(2.3, d)
out <- hitrun(alpha, nbatch = 30, blen = 250,
    a1 = a1, b1 = b1, a2 = a2, b2 = b2)
# prior means
round(colMeans(out$batch), 3)
# Monte Carlo standard errors
round(apply(out$batch, 2, sd) / sqrt(out$nbatch), 3)
}
\keyword{misc}
