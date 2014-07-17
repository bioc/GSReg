/* Registration of C routines */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void Nij( double *vect, int *np,int *mp, double *N);
void kendalltaudist( double *vect, int *np,int *mp, double *dist);
void vect2compC( double *vect, int *np,int *mp, double *comp);



static const R_CMethodDef R_CDef[] = {
  {"Nij", (DL_FUNC)&Nij,4},
  {"kendalltaudist", (DL_FUNC)&kendalltaudist,4},
  {"vect2compC", (DL_FUNC)&vect2compC,4},
  {NULL, NULL, 0},
};

void R_init_GSReg(DllInfo *info)
{
  R_registerRoutines(info,R_CDef,NULL,NULL,NULL);
}
