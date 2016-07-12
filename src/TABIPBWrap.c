/*
 * C routine to interface tabipb with apbs
 *
 * C version authored by:
 * Jiahui Chen, Southern Methodist University, Dallas, TX
 *
 * Additional modifications and updates by:
 * Leighton Wilson, University of Michigan, Ann Arbor, MI
 *
 * Based on package originally written in FORTRAN by:
 * Weihua Geng, Southern Methodist University, Dallas, TX
 * Robery Krasny, University of Michigan, Ann Arbor, MI
 *
 * Rebuild the architecture of wrapper at 6/29/2016
 */

#include "TABIPBstruct.h"

#ifdef TABIPB_APBS
  #include "generic/valist.h"
#endif

int apbs2tabipb_(TABIPBparm* parm, TABIPBvars* vars){

  FILE *wfp;
  char fname_tp[256];
  Vatom *atom;
  int i, ierr;
  extern int tabipb();

  //sprintf(fname_tp, "%s%s.xyzr",parm->fpath, parm->fname);
  sprintf(fname_tp, "molecule.xyzr");
  wfp=fopen(fname_tp,"w");
  for (i = 0; i < parm->number_of_lines; i++) {
    fprintf(wfp, "%f %f %f %f\n", vars->chrpos[3*i], vars->chrpos[3*i + 1],
            vars->chrpos[3*i + 2], vars->atmrad[i]);
  }
  fclose(wfp);

  ierr=tabipb(parm,vars);

  ierr=output_print(vars);

  return 0;
}
