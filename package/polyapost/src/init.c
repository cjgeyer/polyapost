
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "polyapost.h"

static R_NativePrimitiveArgType cwpolya_types[4] =
    {REALSXP, REALSXP, INTSXP, INTSXP};

static R_NativePrimitiveArgType means_types[9] =
    {REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP,
    REALSXP, REALSXP};

static R_NativePrimitiveArgType probvect1_types[8] =
    {REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP};

static R_CMethodDef cMethods[] = {
    {"cwpolya", (DL_FUNC) &cwpolya, 4, cwpolya_types},
    {"means", (DL_FUNC) &means, 9, means_types},
    {"probvect1", (DL_FUNC) &probvect1, 8, probvect1_types},
    {NULL, NULL, 0, NULL}
};

static R_CallMethodDef callMethods[]  = {
    {"hitrun", (DL_FUNC) &hitrun, 11},
    {NULL, NULL, 0}
};

void attribute_visible R_init_polyapost(DllInfo *info)
{
    R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}

