
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include "polyapost.h"

static R_NativePrimitiveArgType cwpolya_types[4] =
    {REALSXP, REALSXP, INTSXP, INTSXP};

static R_NativeArgStyle cwpolya_styles[4] =
    {R_ARG_IN_OUT, R_ARG_IN, R_ARG_IN, R_ARG_IN};

static R_NativePrimitiveArgType means_types[9] =
    {REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP,
    REALSXP, REALSXP};

static R_NativeArgStyle means_styles[9] =
    {R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN,
    R_ARG_IN, R_ARG_IN, R_ARG_OUT};

static R_NativePrimitiveArgType probvect1_types[8] =
    {REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP};

static R_NativeArgStyle probvect1_styles[8] =
    {R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN, R_ARG_IN,
    R_ARG_IN, R_ARG_OUT};

static R_CMethodDef cMethods[] = {
    {"cwpolya", (DL_FUNC) &cwpolya, 4, cwpolya_types, cwpolya_styles},
    {"means", (DL_FUNC) &means, 9, means_types, means_styles},
    {"probvect1", (DL_FUNC) &probvect1, 8, probvect1_types, probvect1_styles},
    {NULL, NULL, 0, NULL, NULL}
};

static R_CallMethodDef callMethods[]  = {
    {"hitrun", (DL_FUNC) &hitrun, 12},
    {NULL, NULL, 0}
};

void R_init_polyapost(DllInfo *info)
{
    R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
}

