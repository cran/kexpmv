#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void mydgexpv_(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mydmexpv_(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void wrapalldgexpv_(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void wrapalldmexpv_(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void wrapdgpadm_(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"mydgexpv_",      (DL_FUNC) &mydgexpv_,      19},
    {"mydmexpv_",      (DL_FUNC) &mydmexpv_,      19},
    {"wrapalldgexpv_", (DL_FUNC) &wrapalldgexpv_, 22},
    {"wrapalldmexpv_", (DL_FUNC) &wrapalldmexpv_, 22},
    {"wrapdgpadm_",    (DL_FUNC) &wrapdgpadm_,    11},
    {NULL, NULL, 0}
};

void R_init_kexpmv(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
