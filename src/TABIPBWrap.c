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

int apbs2tabipb_(TABIPBparm* parm, Valist* molecules){

  FILE *wfp;
  char fname_tp[256];
  Vatom *atom;
  double *t_chrpos, *t_atmchr, *t_atmrad;
  int i, ierr;
  extern int tabipb();

  if ((t_chrpos = (double *) malloc(3 * parm->number_of_lines * sizeof(double))) == NULL) {
          printf("Error in allocating t_chrpos!\n");
  }
  if ((t_atmchr = (double *) malloc(parm->number_of_lines * sizeof(double))) == NULL) {
          printf("Error in allocating t_atmchr!\n");
  }
  if ((t_atmrad = (double *) malloc(parm->number_of_lines * sizeof(double))) == NULL) {
          printf("Error in allocating t_atmrad!\n");
  }

  for (i = 0; i < parm->number_of_lines; i++){
    atom = Valist_getAtom(molecules, i);
    t_chrpos[3*i] = Vatom_getPosition(atom)[0];
    t_chrpos[3*i + 1] = Vatom_getPosition(atom)[1];
    t_chrpos[3*i + 2] = Vatom_getPosition(atom)[2];
    t_atmchr[i] = Vatom_getCharge(atom);
    t_atmrad[i] = Vatom_getRadius(atom);
  }

  sprintf(fname_tp, "%s%s.xyzr",parm->fpath, parm->fname);
  wfp=fopen(fname_tp,"w");
  for (i = 0; i < parm->number_of_lines; i++) {
    fprintf(wfp, "%f %f %f %f\n",t_chrpos[3*i], t_chrpos[3*i + 1], t_chrpos[3*i + 2],
                                 t_atmrad[i]);
  }
  fclose(wfp);

  ierr=tabipb(parm,t_chrpos,t_atmchr,t_atmrad);

  free(t_chrpos);
  free(t_atmchr);
  free(t_atmrad);

  return 0;
}
