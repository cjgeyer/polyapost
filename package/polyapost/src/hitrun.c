
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <string.h>

void propose(double *x, double *proposal, double *a, double *b, int d, int n,
    double *z, double *smax_out, double *smin_out, double *u_out);

static int my_dim_oc = 0;
static int my_dim_nc = 0;
static double *my_alpha = 0;
static double *my_origin = 0;
static double *my_basis = 0;
static double *my_buffer = 0;

static void logh_setup(double *alpha, double *origin, double *basis,
    int dim_oc, int dim_nc)
{
    my_dim_oc = dim_oc;
    my_dim_nc = dim_nc;
    my_alpha = (double *) R_alloc(dim_oc, sizeof(double));
    my_origin = (double *) R_alloc(dim_oc, sizeof(double));
    my_basis = (double *) R_alloc(dim_oc * dim_nc, sizeof(double));
    my_buffer = (double *) R_alloc(dim_oc, sizeof(double));

    memcpy(my_alpha, alpha, dim_oc * sizeof(double));
    memcpy(my_origin, origin, dim_oc * sizeof(double));
    memcpy(my_basis, basis, dim_oc * dim_nc * sizeof(double));
}

static double logh(double *state)
{
    double one = 1.0;
    int ione = 1;
    // my_buffer := my_origin + my_basis * state
    memcpy(my_buffer, my_origin, my_dim_oc * sizeof(double));
    F77_CALL(dgemv)("n", &my_dim_oc, &my_dim_nc, &one,
        my_basis, &my_dim_oc, state, &ione, &one, my_buffer, &ione);

    // return sum((alpha - 1) * my_buffer)
    double result = 0.0;
    for (int i = 0; i < my_dim_oc; i++) {
        if (my_buffer[i] <= 0.0)
            return R_NegInf;
        result += (my_alpha[i] - 1.0) * log(my_buffer[i]);
    }
    return result;
}

static int nrow_my_out_mat = 0;
static int ncol_my_out_mat = 0;
static double *my_out_mat = 0;
static double *my_out_vec = 0;

static void out_setup(double *origin, double *basis, double *outmat,
    int dim_oc, int dim_nc, int dim_out, int has_outmat)
{
    if (has_outmat) {
        my_out_mat = (double *) R_alloc(dim_out * dim_nc, sizeof(double));
        my_out_vec = (double *) R_alloc(dim_out, sizeof(double));
        nrow_my_out_mat = dim_out;
        ncol_my_out_mat = dim_nc;

        // my_out_mat := out_mat * basis
        double one = 1.0;
        double zero = 0.0;
        int ione = 1;
        F77_CALL(dgemm)("n", "n", &dim_out, &dim_nc, &dim_oc, &one, outmat,
            &dim_out, basis, &dim_oc, &zero, my_out_mat, &dim_out);
        // my_out_vec := out_mat * origin
        F77_CALL(dgemv)("n", &dim_out, &dim_nc, &one, outmat, &dim_out,
            origin, &ione, &zero, my_out_vec, &ione);
    } else {
        my_out_mat = (double *) R_alloc(dim_oc * dim_nc, sizeof(double));
        my_out_vec = (double *) R_alloc(dim_oc, sizeof(double));
        nrow_my_out_mat = dim_oc;
        ncol_my_out_mat = dim_nc;
        memcpy(my_out_mat, basis, dim_oc * dim_nc * sizeof(double));
        memcpy(my_out_vec, origin, dim_oc * sizeof(double));
    }
}

static void outfun(double *state, double *buffer)
{
    double one = 1.0;
    int ione = 1;
    // buffer := my_out_vec + my_out_mat * state
    memcpy(buffer, my_out_vec, nrow_my_out_mat * sizeof(double));
    F77_CALL(dgemv)("n", &nrow_my_out_mat, &ncol_my_out_mat, &one,
        my_out_mat, &nrow_my_out_mat, state, &ione, &one, buffer, &ione);
}

void check_finite(double *x, int length, char *name)
{
    for (int i = 0; i < length; i++) {
        if (! R_FINITE(x[i]))
            error("argument \"%s\" must have all components finite", name);
    }
}

void check_positive(double *x, int length, char *name)
{
    for (int i = 0; i < length; i++) {
        if (x[i] <= 0.0)
            error("argument \"%s\" must have all components positive", name);
    }
}

SEXP hitrun(SEXP alpha, SEXP initial, SEXP nbatch, SEXP blen, SEXP nspac,
    SEXP origin, SEXP basis, SEXP amat, SEXP bvec, SEXP outmat, SEXP debug)
{
    if (! isReal(alpha))
        error("argument \"alpha\" must be type double");
    if (! isReal(initial))
        error("argument \"initial\" must be type double");
    if (! isInteger(nbatch))
        error("argument \"nbatch\" must be type integer");
    if (! isInteger(blen))
        error("argument \"blen\" must be type integer");
    if (! isInteger(nspac))
        error("argument \"nspac\" must be type integer");
    if (! isReal(origin))
        error("argument \"origin\" must be type double");
    if (! isReal(basis))
        error("argument \"basis\" must be type double");
    if (! isReal(amat))
        error("argument \"amat\" must be type double");
    if (! isReal(bvec))
        error("argument \"bvec\" must be type double");
    if (! (isNull(outmat) | isReal(outmat)))
        error("argument \"outmat\" must be type double or NULL");
    if (! isLogical(debug))
        error("argument \"debug\" must be logical");

#ifdef BLATHER
    Rprintf("finished checking types\n");
#endif /* BLATHER */

    if (! isMatrix(basis))
        error("argument \"basis\" must be matrix");
    if (! isMatrix(amat))
        error("argument \"amat\" must be matrix");
    if (! (isNull(outmat) | isMatrix(outmat)))
        error("argument \"outmat\" must be matrix or NULL");

#ifdef BLATHER
    Rprintf("finished checking matrices are matrices\n");
#endif /* BLATHER */

    int dim_oc = LENGTH(alpha);
    int dim_nc = LENGTH(initial);
    int ncons = nrows(amat);
    if (LENGTH(nbatch) != 1)
        error("argument \"nbatch\" must be scalar");
    if (LENGTH(blen) != 1)
        error("argument \"blen\" must be scalar");
    if (LENGTH(nspac) != 1)
        error("argument \"nspac\" must be scalar");
    if (LENGTH(origin) != dim_oc)
        error("length(origin) != length(alpha)");
    if (nrows(basis) != dim_oc)
        error("nrow(basis) != length(alpha)");
    if (ncols(basis) != dim_nc)
        error("ncol(basis) != length(initial)");
    if (ncols(amat) != dim_nc)
        error("ncol(amat) != length(initial)");
    if (LENGTH(bvec) != ncons)
        error("length(bvec) != nrow(amat)");
    if (LENGTH(debug) != 1)
        error("argument \"debug\" must be scalar");

    int dim_out = dim_oc;
    if (! isNull(outmat)) {
        dim_out = nrows(outmat);
        if (ncols(outmat) != dim_nc)
            error("ncol(outmat) != length(initial)");
    }

#ifdef BLATHER
    Rprintf("finished checking length of vectors and dimensions of matrices\n");
#endif /* BLATHER */

    int int_nbatch = INTEGER(nbatch)[0];
    int int_blen = INTEGER(blen)[0];
    int int_nspac = INTEGER(nspac)[0];
    int int_debug = LOGICAL(debug)[0];
    double *dbl_star_alpha = REAL(alpha);
    double *dbl_star_initial = REAL(initial);
    double *dbl_star_origin = REAL(origin);
    double *dbl_star_basis = REAL(basis);
    double *dbl_star_amat = REAL(amat);
    double *dbl_star_bvec = REAL(bvec);
    int has_outmat = isMatrix(outmat);
    double *dbl_star_outmat = 0;
    if (has_outmat)
        dbl_star_outmat = REAL(outmat);

    if (int_nbatch <= 0)
        error("argument \"nbatch\" must be positive");
    if (int_blen <= 0)
        error("argument \"blen\" must be positive");
    if (int_nspac <= 0)
        error("argument \"nspac\" must be positive");
    check_finite(dbl_star_alpha, dim_oc, "alpha");
    check_positive(dbl_star_alpha, dim_oc, "alpha");
    check_finite(dbl_star_initial, dim_nc, "initial");
    check_finite(dbl_star_origin, dim_oc, "origin");
    check_finite(dbl_star_basis, dim_oc * dim_nc, "basis");
    check_finite(dbl_star_amat, ncons * dim_nc, "amat");
    check_finite(dbl_star_bvec, ncons * dim_nc, "bvec");
    if (has_outmat)
        check_finite(dbl_star_outmat, ncons * dim_nc, "outmat");

#ifdef BLATHER
    Rprintf("finished checking finiteness and positivity\n");
#endif /* BLATHER */

    double *state = (double *) R_alloc(dim_nc, sizeof(double));
    double *proposal = (double *) R_alloc(dim_nc, sizeof(double));
    double *batch_buffer = (double *) R_alloc(dim_out, sizeof(double));
    double *out_buffer = (double *) R_alloc(dim_out, sizeof(double));

    memcpy(state, dbl_star_initial, dim_nc * sizeof(double));
    logh_setup(dbl_star_alpha, dbl_star_origin, dbl_star_basis, dim_oc, dim_nc);
    double current_log_dens = logh(state);

    out_setup(dbl_star_origin, dbl_star_basis, dbl_star_outmat, dim_oc, dim_nc,
        dim_out, has_outmat);

    SEXP result, resultnames, path, save_initial, save_final;

    if (! int_debug) {
        PROTECT(result = allocVector(VECSXP, 3));
        PROTECT(resultnames = allocVector(STRSXP, 3));
    } else {
        PROTECT(result = allocVector(VECSXP, 11));
        PROTECT(resultnames = allocVector(STRSXP, 11));
    }
    PROTECT(path = allocMatrix(REALSXP, dim_out, int_nbatch));
    SET_VECTOR_ELT(result, 0, path);
    PROTECT(save_initial = duplicate(initial));
    SET_VECTOR_ELT(result, 1, save_initial);
    UNPROTECT(2);
    SET_STRING_ELT(resultnames, 0, mkChar("batch"));
    SET_STRING_ELT(resultnames, 1, mkChar("initial"));
    SET_STRING_ELT(resultnames, 2, mkChar("final"));
    if (int_debug) {
        SEXP spath, ppath, zpath, u1path, u2path, s1path, s2path, gpath;
        int nn = int_nbatch * int_blen * int_nspac;
        PROTECT(spath = allocMatrix(REALSXP, dim_nc, nn));
        SET_VECTOR_ELT(result, 3, spath);
        PROTECT(ppath = allocMatrix(REALSXP, dim_nc, nn));
        SET_VECTOR_ELT(result, 4, ppath);
        PROTECT(zpath = allocMatrix(REALSXP, dim_nc, nn));
        SET_VECTOR_ELT(result, 5, zpath);
        PROTECT(u1path = allocVector(REALSXP, nn));
        SET_VECTOR_ELT(result, 6, u1path);
        PROTECT(u2path = allocVector(REALSXP, nn));
        SET_VECTOR_ELT(result, 7, u2path);
        PROTECT(s1path = allocVector(REALSXP, nn));
        SET_VECTOR_ELT(result, 8, s1path);
        PROTECT(s2path = allocVector(REALSXP, nn));
        SET_VECTOR_ELT(result, 9, s2path);
        PROTECT(gpath = allocVector(REALSXP, nn));
        SET_VECTOR_ELT(result, 10, gpath);
        UNPROTECT(8);
        SET_STRING_ELT(resultnames, 3, mkChar("current"));
        SET_STRING_ELT(resultnames, 4, mkChar("proposal"));
        SET_STRING_ELT(resultnames, 5, mkChar("z"));
        SET_STRING_ELT(resultnames, 6, mkChar("u1"));
        SET_STRING_ELT(resultnames, 7, mkChar("u2"));
        SET_STRING_ELT(resultnames, 8, mkChar("s1"));
        SET_STRING_ELT(resultnames, 9, mkChar("s2"));
        SET_STRING_ELT(resultnames, 10, mkChar("log.green"));
    }
    namesgets(result, resultnames);
    UNPROTECT(1);

    // REVISED DOWN TO HERE

    GetRNGstate();

    if (current_log_dens == R_NegInf)
        error("log unnormalized density -Inf at initial state");

    for (int ibatch = 0, k = 0; ibatch < int_nbatch; ibatch++) {

        int j;

        for (int i = 0; i < dim_out; i++)
            batch_buffer[i] = 0.0;

        for (int jbatch = 0; jbatch < int_blen; jbatch++) {

            double proposal_log_dens;

            for (int ispac = 0; ispac < int_nspac; ispac++) {

                /* Note: should never happen! */
                if (current_log_dens == R_NegInf)
                    error("log density -Inf at current state");

                double u1 = R_NaReal;
                double u2 = R_NaReal;
                double smax = R_NaReal;
                double smin = R_NaReal;
                double z[dim_nc];

                propose(state, proposal, dbl_star_amat, dbl_star_bvec,
                    dim_nc, ncons, z, &smax, &smin, &u1);

                proposal_log_dens = logh(proposal);

                int accept = FALSE;
                if (proposal_log_dens != R_NegInf) {
                    if (proposal_log_dens >= current_log_dens) {
                        accept = TRUE;
                    } else {
                        double green = exp(proposal_log_dens
                            - current_log_dens);
                        u2 = unif_rand();
                        accept = u2 < green;
                    }
                }

                if (int_debug) {
                    int l = ispac + int_nspac * (jbatch + int_blen * ibatch);
                    int lbase = l * dim_nc;
                    SEXP spath = VECTOR_ELT(result, 3);
                    SEXP ppath = VECTOR_ELT(result, 4);
                    SEXP zpath = VECTOR_ELT(result, 5);
                    SEXP u1path = VECTOR_ELT(result, 6);
                    SEXP u2path = VECTOR_ELT(result, 7);
                    SEXP s1path = VECTOR_ELT(result, 8);
                    SEXP s2path = VECTOR_ELT(result, 9);
                    SEXP gpath = VECTOR_ELT(result, 10);
                    for (int lj = 0; lj < dim_nc; lj++) {
                        REAL(spath)[lbase + lj] = state[lj];
                        REAL(ppath)[lbase + lj] = proposal[lj];
                        REAL(zpath)[lbase + lj] = z[lj];
                    }
                    REAL(u1path)[l] = u1;
                    REAL(u2path)[l] = u2;
                    REAL(s1path)[l] = smin;
                    REAL(s2path)[l] = smax;
                    REAL(gpath)[l] = proposal_log_dens - current_log_dens;
                }

                if (accept) {
                    memcpy(proposal, state, dim_nc * sizeof(double));
                    current_log_dens = proposal_log_dens;
                }
            } /* end of inner loop (one iteration) */

            outfun(state, out_buffer);
            for (j = 0; j < dim_out; j++)
                batch_buffer[j] += out_buffer[j];

        } /* end of middle loop (one batch) */

        for (j = 0; j < dim_out; j++, k++)
            REAL(path)[k] = batch_buffer[j] / int_blen;

    } /* end of outer loop */

    PutRNGstate();

    PROTECT(save_final = allocVector(REALSXP, dim_nc));
    memcpy(REAL(save_final), state, dim_nc * sizeof(double));
    SET_VECTOR_ELT(result, 2, save_final);

    UNPROTECT(5);
    return result;
}

void propose(double *x, double *proposal, double *a, double *b, int d, int n,
    double *z, double *smax_out, double *smin_out, double *u_out)
{
    for (int i = 0; i < d; i++) {
        z[i] = norm_rand();
    }

    double smax = R_PosInf;
    double smin = R_NegInf;

    for (int i = 0; i < n; i++) {

        double ax = 0.0;
        double az = 0.0;
        for (int j = 0; j < d; j++) {
            ax += a[i + j * n] * x[j];
            az += a[i + j * n] * z[j];
        }
        double bound = (b[i] - ax) / az;
        if (az > 0 && bound < smax)
                smax = bound;
        if (az < 0 && bound > smin)
                smin = bound;
    }

    double u = unif_rand();

    for (int i = 0; i < d; i++)
        proposal[i] = x[i] + (u * smin + (1.0 - u) * smax) * z[i];

    *smax_out = smax;
    *smin_out = smin;
    *u_out = u;
}

