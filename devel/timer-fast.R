
 library(polyapost, lib.loc = "../package/polyapost.Rcheck")

 d <- 601
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
 #     constraints works for both (except this takes forever when
 #     dimension is big) so we actually do case split
 # mean constraint
 a2 <- rbind(x)
 b2 <- (d + 1) / 2
 # median constraint(s)
 if (d %% 2 == 0) {
     a1 <- NULL
     b1 <- NULL
     a2 <- rbind(a2, as.numeric(2 * x < d + 1))
     b2 <- c(b2, 1 / 2)
 } else {
     a1 <- rbind(- as.numeric(2 * x <= d + 1), - as.numeric(2 * x >= d + 1))
     b1 <- c(- 1 / 2, - 1 / 2)
 }
 # median absolute deviation constraint
 r <- floor((d - 1) / 4)
 a1 <- rbind(a1, - as.numeric(x <= (d + 1) / 2 - r) -
    as.numeric(x >= (d + 1) / 2 + r))
 b1 <- c(b1, - 1 / 2)

 # simulate prior, which is Dirichlet(alpha)
 # posterior would be another Dirichlet with n + alpha - 1,
 #    where n is count of IID data for each value
 alpha <- rep(2.3, d)
 out <- hitrun(alpha, nbatch = 30, blen = 250,
     a1 = a1, b1 = b1, a2 = a2, b2 = b2)
 out$time
 out$split.time

