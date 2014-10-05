
# a1 and b1 define inequality constraints (a1 %*% x <= b1)
# a2 and b2 define inequality constraints (a2 %*% x == b2)
# as in the makeH function in the rcdd package
# otherwise arguments are like metrop in the mcmc package
# we add to a1 and b1 the unit simplex constraints (all(x >= 0))
# we add to a2 and b2 the unit simplex constraints (sum(x) == 1)
# a1 and b1 can be missing
# a2 and b2 can be missing

hitrun <- function(alpha, ...)
    UseMethod("hitrun")

hitrun.hitrun <- function(alpha, nbatch, blen, nspac, outmat, debug, ...)
{
    if (missing(nbatch)) nbatch <- alpha$nbatch
    if (missing(blen)) blen <- alpha$blen
    if (missing(nspac)) nspac <- alpha$nspac
    if (missing(outmat)) outmat <- alpha$outmat
    if (missing(debug)) debug <- alpha$debug

    assign(".Random.seed", alpha$final.seed, .GlobalEnv)
    basis <- alpha$basis
    origin <- alpha$origin
    amat <- alpha$amat
    bvec <- alpha$bvec

    out <- hitrunHelper(alpha$alpha, alpha$final, nbatch, blen, nspac,
        origin, basis, amat, bvec, outmat, debug)

    class(out) <- "hitrun"
    return(out)
}

hitrun.default <- function(alpha, a1 = NULL, b1 = NULL, a2 = NULL, b2 = NULL, 
     nbatch = 1, blen = 1, nspac = 1, outmat = NULL, debug = FALSE, ...)
{
    if (! exists(".Random.seed")) runif(1)
    saveseed <- .Random.seed

    stopifnot(is.numeric(alpha))
    stopifnot(is.finite(alpha))
    stopifnot(is.vector(alpha))
    stopifnot(length(alpha) >= 2)
    stopifnot(alpha > 0)
    stopifnot(is.null(a1) == is.null(b1))
    if (! is.null(a1)) {
        stopifnot(is.numeric(a1))
        stopifnot(is.finite(a1))
        stopifnot(is.matrix(a1))
        stopifnot(is.numeric(b1))
        stopifnot(is.finite(b1))
        stopifnot(is.vector(b1))
        stopifnot(nrow(a1) == length(b1))
        stopifnot(ncol(a1) == length(alpha))
    }
    stopifnot(is.null(a2) == is.null(b2))
    if (! is.null(a2)) {
        stopifnot(is.numeric(a2))
        stopifnot(is.finite(a2))
        stopifnot(is.matrix(a2))
        stopifnot(is.numeric(b2))
        stopifnot(is.finite(b2))
        stopifnot(is.vector(b2))
        stopifnot(nrow(a2) == length(b2))
        stopifnot(ncol(a2) == length(alpha))
    }

    if (! is.null(outmat)) {
        stopifnot(is.numeric(outmat))
        stopifnot(is.finite(outmat))
        stopifnot(is.matrix(outmat))
        stopifnot(ncol(outmat) == length(alpha))
    }

    stopifnot(is.numeric(nbatch))
    stopifnot(nbatch == as.integer(nbatch))
    stopifnot(length(nbatch) == 1)
    stopifnot(nbatch >= 1)
    stopifnot(is.numeric(blen))
    stopifnot(blen == as.integer(blen))
    stopifnot(length(blen) == 1)
    stopifnot(blen >= 1)
    stopifnot(is.numeric(nspac))
    stopifnot(nspac == as.integer(nspac))
    stopifnot(length(nspac) == 1)
    stopifnot(nspac >= 1)
    stopifnot(is.logical(debug))
    stopifnot(length(debug) == 1)

    d <- length(alpha)
    # add unit simpled constraints
    a1 <- rbind(a1, - diag(d))
    b1 <- c(b1, rep(0, d))
    a2 <- rbind(a2, rep(1, d))
    b2 <- c(b2, 1)

    hrep1 <- makeH(a1, b1, a2, b2)
    hrep1 <- d2q(hrep1)
    hrep2 <- redundant(hrep1)$output
    hrep3 <- hrep2[hrep2[ , 1] == "1", , drop = FALSE]
    hrep4 <- hrep2[hrep2[ , 1] == "0", , drop = FALSE]

    vrep3 <- scdd(hrep3, representation = "H")$output
    is.line <- vrep3[ , 1] == "1" & vrep3[ , 2] == "0"
    is.point <- vrep3[ , 1] == "0" & vrep3[ , 2] == "1"
    if (! all(is.point | is.line))
        stop("unexpected V-representation of affine hull of constraint set")
    if (sum(is.point) != 1)
        stop("unexpected V-representation of affine hull of constraint set")
    if (sum(is.line) != ncol(a1) - nrow(hrep3))
        stop("unexpected V-representation of affine hull of constraint set")

    foo <- vrep3[ , - c(1, 2), drop = FALSE]
    origin <- foo[is.point, ]
    basis <- foo[is.line, , drop = FALSE]
    basis <- t(basis)

    # at this point
    #     fred <- function(x) origin + basis %*% x
    # maps from new coordinates (NC) onto the affine hull
    #     of the constraint (a convex polytope) in original coordinated (OC)

    amat <- qneg(hrep4[ , - c(1, 2), drop = FALSE])
    bvec <- hrep4[ , 2]
    bvec <- qmq(bvec, qmatmult(amat, cbind(origin)))
    amat <- qmatmult(amat, basis)

    # at this point
    #     sally <- function(x) all(amat %*% x <= bvec)
    # is the indicator function of a convex polytope in NC that is
    # mapped one-to-one onto the constraint set in OC by the function
    # fred defined in the previous comment

    hrep5 <- cbind("0", bvec, qneg(amat))
    vrep5 <- scdd(hrep5, representation = "H")$output
    if (nrow(vrep5) == 0)
        stop("empty constraint set")
    is.point <- vrep5[ , 1] == "0" & vrep5[ , 2] == "1"
    if (! all(is.point))
        stop("unbounded constraint set (can't happen)")
    v <- vrep5[ , - c(1, 2), drop = FALSE]
    rip <- apply(v, 2, qsum)
    rip <- qdq(rip, rep(as.character(nrow(v)), length(rip)))

    # rip is relative interior point of constraint set in NC

    origin <- q2d(origin)
    basis <- q2d(basis)
    bvec <- q2d(bvec)
    amat <- q2d(amat)
    rip <- q2d(rip)

    out <- hitrunHelper(alpha, rip, nbatch, blen, nspac,
        origin, basis, amat, bvec, outmat, debug)

    if (! missing(a1)) {
        out$a1 <- a1
        out$b1 <- b1
    }
    if (! missing(a2)) {
        out$a2 <- a2
        out$b2 <- b2
    }
    class(out) <- "hitrun"
    return(out)
}

hitrunHelper <- function(alpha, initial, nbatch, blen, nspac,
        origin, basis, amat, bvec, outmat, debug)
{
    if (! exists(".Random.seed")) runif(1)
    saveseed <- .Random.seed

    stopifnot(is.numeric(alpha))
    stopifnot(is.finite(alpha))
    stopifnot(is.vector(alpha))
    stopifnot(length(alpha) >= 2)
    stopifnot(alpha > 0)

    stopifnot(is.numeric(initial))
    stopifnot(is.finite(initial))

    stopifnot(is.numeric(nbatch))
    stopifnot(length(nbatch) == 1)
    stopifnot(nbatch == round(nbatch))
    stopifnot(nbatch >= 1)
    stopifnot(is.numeric(blen))
    stopifnot(length(blen) == 1)
    stopifnot(blen == round(blen))
    stopifnot(blen >= 1)
    stopifnot(is.numeric(nspac))
    stopifnot(length(nspac) == 1)
    stopifnot(nspac == round(nspac))
    stopifnot(nspac >= 1)

    stopifnot(is.numeric(basis))
    stopifnot(is.finite(basis))
    stopifnot(is.matrix(basis))
    stopifnot(is.numeric(origin))
    stopifnot(is.finite(origin))
    stopifnot(is.vector(origin))
    stopifnot(nrow(basis) == length(alpha))
    stopifnot(ncol(basis) == length(initial))
    stopifnot(length(origin) == length(alpha))

    stopifnot(is.numeric(amat))
    stopifnot(is.finite(amat))
    stopifnot(is.matrix(amat))
    stopifnot(is.numeric(bvec))
    stopifnot(is.finite(bvec))
    stopifnot(is.vector(bvec))
    stopifnot(nrow(amat) == length(bvec))
    stopifnot(ncol(amat) == length(initial))

    if (! is.null(outmat)) {
        stopifnot(is.numeric(outmat))
        stopifnot(is.finite(outmat))
        stopifnot(is.matrix(outmat))
        stopifnot(ncol(outmat) == length(alpha))
    }

    stopifnot(is.logical(debug))
    stopifnot(length(debug) == 1)

    storage.mode(basis) <- "double"
    storage.mode(amat) <- "double"
    if (is.matrix(outmat))
        storage.mode(outmat) <- "double"
    out.time <- system.time(
    out <- .Call(C_hitrun, as.double(alpha), as.double(initial),
        as.integer(nbatch), as.integer(blen), as.integer(nspac),
        as.double(origin), basis, amat, as.double(bvec), outmat,
        as.logical(debug))
    )

    out$initial.seed <- saveseed
    out$final.seed <- .Random.seed
    out$time <- out.time
    out$alpha <- alpha
    out$nbatch <- nbatch
    out$blen <- blen
    out$nspac <- nspac
    out$origin <- origin
    out$basis <- basis
    out$amat <- amat
    out$bvec <- bvec
    out$outmat <- outmat
    out$batch <- t(out$batch)
    out$debug <- debug
    if (debug) {
        out$current <- t(out$current)
        out$proposal <- t(out$proposal)
        out$z <- t(out$z)
    }
    class(out) <- "hitrun"
    return(out)
}

