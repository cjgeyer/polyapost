
 library(polyapost, lib.loc = "../package/polyapost.Rcheck")

 d <- 11
 x <- 1:d
 # equality constraints
 #     mean equal to (d + 1) / 2
 #     median equal to (d + 1) / 2
 # inequality constraints
 #     median absolute deviation to be at least floor((d - 1) / 4)
 # note: median constraint is tricky.  When d is even it requires
 #         Pr(X > (d + 1) / 2) = 1 / 2
 #     but when d is odd it requires
 #         Pr(X >= (d + 1) / 2) >= 1 / 2
 #     and
 #         Pr(X <= (d + 1) / 2) >= 1 / 2
 #     so is really two inequality constraints rather than
 #     one equality constraint
 #     However, we do not need to case split on whether d is odd
 #     or even because the characterization with two inequality
 #     constraints works for both
 a1 <- rbind(as.numeric(2 * x <= d + 1), as.numeric(2 * x >= d + 1))
 b1 <- c(1 / 2, 1 / 2)
 r <- floor((d - 1) / 4)
 a1 <- rbind(a1, - as.numeric(x <= (d + 1) / 2 - r) -
    as.numeric(x >= (d + 1) / 2 + r))
 b1 <- c(b1, - 1 / 2)
 a2 <- rbind(x)
 b2 <- (d + 1) / 2
 # simulate prior, which is Dirichlet(alpha)
 # posterior would be another Dirichlet with n + alpha - 1,
 #    where n is count of IID data for each value
 alpha <- rep(2.3, d)
 try(hitrun(alpha, nbatch = 30, blen = 250,
     a1 = a1, b1 = b1, a2 = a2, b2 = b2))
 ##### WTF ?????
 a1
 b1
 a2
 b2
 ## Oh!  We have the wrong constraints.  But still should give sane error.

 a1 <- d2q(a1)
 b1 <- d2q(b1)
 a2 <- d2q(a2)
 b2 <- d2q(b2)
 a1 <- rbind(a1, d2q(- diag(d)))
 b1 <- c(b1, rep("0", d))
 a2 <- rbind(a2, rep("1", d))
 b2 <- c(b2, "1")

 hrep1 <- makeH(a1, b1, a2, b2)
 linout <- linearity(hrep1)
 linout
 ## why is constraint 11 an equality constraint?
 a11 <- qneg(hrep1[11, - c(1, 2)])
 b11 <- hrep1[11, 2]
 ## constraint is a11^T x <= b11
 a11
 b11
 ## Oh.  I see!  The only way we can have Pr(X <= 6) <= 1 / 2 and
 ## Pr(X >= 6) <= 1 / 2 is if Pr(X == 6) == 0

 hrep1[linout, 1] <- "1"
 hrep3 <- hrep1[hrep1[ , 1] == "1", , drop = FALSE]
 hrep4 <- hrep1[hrep1[ , 1] == "0", , drop = FALSE]

 hrep3
 hrep4

 vrep3 <- scdd(hrep3, representation = "H")$output
 vrep3

 is.line <- vrep3[ , 1] == "1" & vrep3[ , 2] == "0"
 is.point <- vrep3[ , 1] == "0" & vrep3[ , 2] == "1"
 if (! all(is.point | is.line))
     stop("unexpected V-representation of affine hull of constraint set")
 if (sum(is.point) != 1)
     stop("unexpected V-representation of affine hull of constraint set")
 if (sum(is.line) == 0)
     stop("constraint set consists of single point or is empty")

 foo <- vrep3[ , - c(1, 2), drop = FALSE]
 origin <- foo[is.point, ]
 basis <- foo[is.line, , drop = FALSE]
 basis <- t(basis)
 origin
 basis

 # does this satisfy original equality constraints and the additional
 # ones we discovered,
 # that is, a2 %*% (origin + basis %*% x) == b2, for all possible x,
 # where we have added the constraints Pr(X < 6) = Pr(X > 6) = 1 / 2
 # and Pr(X = 6) = 0?
 # And that is equivalent to a2 %*% origin == b2 and a2 %*% basis == 0
 a2 <- rbind(a2, as.character(as.numeric(x < 6)),
     as.character(as.numeric(x > 6)), as.character(as.numeric(x == 6)))
 b2 <- c(b2, "1/2", "1/2", "0")
 qmq(b2, qmatmult(a2, cbind(origin)))
 qmatmult(a2, basis)
 # seems OK to here

 amat <- qneg(hrep4[ , - c(1, 2), drop = FALSE])
 bvec <- hrep4[ , 2]
 bvec <- qmq(bvec, qmatmult(amat, cbind(origin)))
 amat <- qmatmult(amat, basis)
 amat
 bvec

 hrep5 <- cbind("0", bvec, qneg(amat), "-1")
 grad5 <- c(rep("0", ncol(amat)), "1")
 lout <- lpcdd(hrep5, grad5, minimize = FALSE)
 if (lout$solution.type != "Optimal" || qsign(lout$optimal.value) <= 0)
     stop("constraint set is empty (constraints are inconsistent)")
 x <- lout$primal.solution
 # rip is relative interior point of constraint set in NC
 rip <- x[- length(x)]
 rip
 # back to original coordinates
 qpq(origin, qmatmult(basis, cbind(rip)))
 # that looks OK too

 origin <- q2d(origin)
 basis <- q2d(basis)
 bvec <- q2d(bvec)
 amat <- q2d(amat)
 rip <- q2d(rip)

 nbatch <- 30
 blen <- 250
 nspac <- 1
 outmat <- NULL
 debug <- FALSE
 try(polyapost:::hitrunHelper(alpha, rip, nbatch, blen, nspac,
        origin, basis, amat, bvec, outmat, debug))

 # so do in computer inexact arithmetic
 xinit <- origin + basis %*% rip
 xinit
 # And that's right.  We do not allow components of x to be non-positive.

 # can we give more informative error message
 hrep3
 a3 <- hrep3[ , - c(1, 2), drop = FALSE]
 b3 <- hrep3[ , 2]
 solo <- apply(a3, 1, function(x) sum(x != "0") == 1)
 solo
 a3[solo, ]
 b3[solo]
 any(b3[solo] == "0")
 # that's it
 molo <- a3[solo, , drop = FALSE]
 nolo <- apply(molo, 1, function(x) seq(along = x)[x != "0"])
 colo <- paste(nolo, collapse = ", ")
 paste("variable(s)", colo, "constrained to be exactly zero")

 # now for the timing

