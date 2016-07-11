/*
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

int main(int argc, char *argv[]){
  /* main reads the input file, writes xyzr file for msms and sets up position,
     raduis and charges */

  FILE *fp, *wfp;
  char c[16];
  char fname_tp[256];
  char c1[10], c2[10], c3[10], c4[10], c5[10];
  double a1, a2, a3, b1, b2, b3;
  double epsp, epsw, bulk_strength, theta, temp;
  int maxparnode,order,ierr,i;
  double *t_chrpos, *t_atmchr, *t_atmrad;

  /* time */
  extern void timer_start();
  extern void timer_end();

  timer_start("TOTAL_TIME");

  TABIPBparm *main_parm;
  main_parm = (TABIPBparm*)calloc(1,sizeof(TABIPBparm));
  TABIPBvars *main_vars;
  main_vars = (TABIPBvars*)calloc(1,sizeof(TABIPBvars));

  extern int tabipb();
/********************************************************/

  fp=fopen("usrdata.in","r");
    ierr=fscanf(fp,"%s %s",c,main_parm->fname);
    ierr=fscanf(fp,"%s %s",c,main_parm->density);
    ierr=fscanf(fp,"%s %s",c,main_parm->probe_radius);
    ierr=fscanf(fp,"%s %lf",c,&epsp);
    main_parm->epsp = epsp;
    ierr=fscanf(fp,"%s %lf",c,&epsw);
    main_parm->epsw = epsw;
    ierr=fscanf(fp,"%s %lf",c,&bulk_strength);
    main_parm->bulk_strength = bulk_strength;
    ierr=fscanf(fp,"%s %d",c,&order);
    main_parm->order = order;
    ierr=fscanf(fp,"%s %d",c,&maxparnode);
    main_parm->maxparnode = maxparnode;
    ierr=fscanf(fp,"%s %lf",c,&theta);
    main_parm->theta = theta;
    ierr=fscanf(fp,"%s %lf",c,&temp);
    main_parm->temp = temp;
  fclose(fp);

  main_parm->temp = 300;

  main_parm->mesh_flag = 0;

/********************************************************/

  sprintf(main_parm->fpath, "");

  sprintf(fname_tp, "%s%s.pqr", main_parm->fpath, main_parm->fname);
  fp = fopen(fname_tp, "r");
  sprintf(fname_tp, "%s%s.xyzr", main_parm->fpath, main_parm->fname);
  wfp=fopen(fname_tp, "w");
  int ch;// main_parm->number_of_lines = 0;

  while(fscanf(fp, "%s %s %s %s %s %lf %lf %lf %lf %lf",
               c1, c2, c3, c4, c5, &a1, &a2, &a3, &b1, &b2) != EOF) {
    fprintf(wfp, "%f %f %f %f\n", a1, a2, a3, b2);
    main_parm->number_of_lines++;
  }

  fclose(fp);
  fclose(wfp);

  sprintf(fname_tp, "%s%s.pqr", main_parm->fpath, main_parm->fname);
  fp = fopen(fname_tp, "r");

  if ((main_vars->chrpos = (double *) malloc(3 * main_parm->number_of_lines * sizeof(double))) == NULL) {
          printf("Error in allocating t_chrpos!\n");
  }
  if ((main_vars->atmchr = (double *) malloc(main_parm->number_of_lines * sizeof(double))) == NULL) {
          printf("Error in allocating t_atmchr!\n");
  }
  if ((main_vars->atmrad = (double *) malloc(main_parm->number_of_lines * sizeof(double))) == NULL) {
          printf("Error in allocating t_atmrad!\n");
  }

  for (i = 0; i < main_parm->number_of_lines; i++) {
    ierr = fscanf(fp, "%s %s %s %s %s %lf %lf %lf %lf %lf",
                  c1,c2,c3,c4,c5,&a1,&a2,&a3,&b1,&b2);
    main_vars->chrpos[3*i] = a1;
    main_vars->chrpos[3*i + 1] = a2;
    main_vars->chrpos[3*i + 2] = a3;
    main_vars->atmchr[i] = b1;
    main_vars->atmrad[i] = b2;
  }

  fclose(fp);
  printf("Finished assembling charge (.pqr) file...\n");

  ierr=tabipb(main_parm, main_vars);

  free(main_parm);
  free(main_vars->atmchr);
  free(main_vars->chrpos);
  free(main_vars->atmrad);
  free(main_vars->vert_ptl); // allocate in output_potential()
  free(main_vars->xvct);
  free(main_vars);

  timer_end();

  return 0;
}
