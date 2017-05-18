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
     nbatch = 1, blen = 1, nspac = 1, outmat = NULL, debug = FALSE,
     stop.if.implied.equalities = FALSE, \dots)

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
  \item{stop.if.implied.equalities}{If \code{TRUE} stop if there
      are any implied equalities.}
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
See also the section on this subject below.
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
\section{GMP Rational Arithmetic}{
The arguments \code{a1}, \code{b1}, \code{a2}, and \code{b2} can and should be
given as GMP (GNU multiple precision) rational values.  This allows the
computational geometry calculations for the constraint set to be done exactly,
without error.  For example, if \code{a1} has elements that have been rounded
to two decimal places one should do
\preformatted{
a1 <- z2q(round(100 * a1), rep(100, length(a1)))
}
and similarly for \code{b1}, \code{a2}, and \code{b2} to make them exact.
For all the conversion functions between ordinary computer numbers and
GMP rational numbers see \link{ConvertGMP}.  For all the functions that
do arithmetic on GMP rational numbers, see \link{ArithmeticGMP}.
}
\section{Warning About Implied Equality Constraints}{
If any constraints supplied as inequality constraints (specified by rows
of \code{a1} and the corresponding components of \code{b1}) actually hold
with equality for all points in the constraint set, this is called an
implied equality constraint.  The program must establish that none of these
exist (which is a fast operation) or, otherwise, find out which constraints
supplied as inequality constraints are actually implied equality constraints,
and this operation is very slow when the state is high dimensional.  One
example with 1000 variables took 3 days of computing time when there were
implied equality constraints in the specification.  The same example takes
9 minutes when the same constraint set is specified in a different way so
that there are no  implied equality constraints.

This issue is not a big deal if there are only in the low hundreds of
variables, because the algorithm to find implied equality constraints is not
that slow.  The same example that takes 3 days of computing time with 1000
variables takes only 15 seconds with 100 variables, 3 and 1/2 minutes with 200
variables, and 23 minutes with 300 variables.  As one can see, this issue
does become a big deal as the number of variables increases.  Thus users
should avoid implied inequality constraints, if possible,
when there are many variables.  Admittedly, there
is no sure way users can identify and eliminate implied equality constraints.
(The sure way to do that is precisely the time consuming step we are trying
to avoid.)  The argument \code{stop.if.implied.equalities} can be used to
quickly test for the presence of implied equalities.
}
\section{Philosophy of MCMC}{
This function follows the philosophy of MCMC used in the CRAN package
\code{mcmc} and the introductory chapter of the
\emph{Handbook of Markov Chain Monte Carlo} (Geyer, 2011).

The \code{hitrun} function automatically does batch means in order to reduce
the size of output and to enable easy calculation of Monte Carlo standard
errors (MCSE), which measure error due to the Monte Carlo sampling (not
error due to statistical sampling --- MCSE gets smaller when you run the
computer longer, but statistical sampling variability only gets smaller
when you get a larger data set).  All of this is explained in the package
vignette for the \code{mcmc} package (\code{vignette("demo", "mcmc")})
and in Section 1.10 of Geyer (2011).

The \code{hitrun} function does not apparently
do \dQuote{burn-in} because this concept does not actually help with MCMC
(Geyer, 2011, Section 1.11.4) but the re-entrant property of the hitrun
function does allow one to do \dQuote{burn-in} if one wants.
Assuming \code{alpha}, \code{a1}, \code{b1}, \code{a2}, and \code{b2}
have been already defined
\preformatted{
out <- hitrun(alpha, a1, b1, a2, b2, nbatch = 1, blen = 1e5)
out <- hitrun(out, nbatch = 100, blen = 1000)
}
throws away a run of 100 thousand iterations before doing another run of
100 thousand iterations that is actually useful for analysis, for example,
\preformatted{
apply(out$batch, 2, mean)
apply(out$batch, 2, sd)
}
gives estimates of posterior means and their MCSE assuming the batch length
(here 1000) was long enough to contain almost all of the signifcant
autocorrelation (see Geyer, 2011, Section 1.10, for more on MCSE).
The re-entrant property of the \code{hitrun} function (the second run starts
where the first one stops) assures that this is really \dQuote{burn-in}.

The re-entrant property allows one to do very long runs without having to
do them in one invocation of the \code{hitrun} function.
\preformatted{
out2 <- hitrun(out)
out3 <- hitrun(out2)
batch <- rbind(out$batch, out2$batch, out3$batch)
}
produces a result as if the first run had been three times as long.
}
\references{
Belisle, C. J. P., Romeijn, H. E. and Smith, R. L. (1993)
Hit-and-run algorithms for generating multivariate distributions.
\emph{Mathematics of Operations Research}, \bold{18}, 255--266.

Chen, M. H. and Schmeiser, B. (1993)
Performance of the Gibbs, hit-and-run, and Metropolis samplers.
\emph{Journal of Computational and Graphical Statistics}, \bold{2}, 251--272.

Geyer, C. J. (2011)
Introduction to MCMC.
In \emph{Handbook of Markov Chain Monte Carlo}. Edited by S. P. Brooks,
A. E. Gelman, G. L. Jones, and X. L. Meng.
Chapman & Hall/CRC, Boca Raton, FL, pp. 3--48.
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
