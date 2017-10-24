
#include "header.h"

#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void GIBBS_gp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GIBBSsp_gp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GIBBSsptp_gp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void GIBBStp_gp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void z_pr_its_gp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void z_pr_its_gp_sp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void z_pr_its_gp_sptp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void z_pr_its_gp_tp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void zlt_fore_gp_its(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void zlt_fore_gp_sp_its(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void zlt_fore_gp_sptp_its(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void zlt_fore_gp_tp_its(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"GIBBS_gp",             (DL_FUNC) &GIBBS_gp,             44},
    {"GIBBSsp_gp",           (DL_FUNC) &GIBBSsp_gp,           50},
    {"GIBBSsptp_gp",         (DL_FUNC) &GIBBSsptp_gp,         64},
    {"GIBBStp_gp",           (DL_FUNC) &GIBBStp_gp,           57},
    {"z_pr_its_gp",          (DL_FUNC) &z_pr_its_gp,          22},
    {"z_pr_its_gp_sp",       (DL_FUNC) &z_pr_its_gp_sp,       28},
    {"z_pr_its_gp_sptp",     (DL_FUNC) &z_pr_its_gp_sptp,     36},
    {"z_pr_its_gp_tp",       (DL_FUNC) &z_pr_its_gp_tp,       30},
    {"zlt_fore_gp_its",      (DL_FUNC) &zlt_fore_gp_its,      22},
    {"zlt_fore_gp_sp_its",   (DL_FUNC) &zlt_fore_gp_sp_its,   26},
    {"zlt_fore_gp_sptp_its", (DL_FUNC) &zlt_fore_gp_sptp_its, 33},
    {"zlt_fore_gp_tp_its",   (DL_FUNC) &zlt_fore_gp_tp_its,   30},
    {NULL, NULL, 0}
};

void R_init_spTDyn(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

//
