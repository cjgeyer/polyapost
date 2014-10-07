
#ifndef POLYAPOST_H
#define POLYAPOST_H

void cwpolya(double *x, double *w, int *nin, int *Nin);

void means(double *p, int *matsize, double *a2, int *nrow, double *b2,
    double *initsol, int *rep, double *ysamp, double *estimate);

void probvect1(double *p, int *matsize, double *a2, int *nrow, double *b2,
    double *initsol, int *length, double *estimate);

#include <Rinternals.h>

SEXP hitrun(SEXP alpha, SEXP initial, SEXP nbatch, SEXP blen, SEXP nspac,
    SEXP origin, SEXP basis, SEXP amat, SEXP bvec, SEXP outmat, SEXP debug,
    SEXP scale);

#endif /* POLYAPOST_H */

