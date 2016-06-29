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

#include <stdlib.h> /* calloc() */
#include <stdio.h> /* FILE */

int main(int argc, char *argv[]){

  FILE *fp, *wfp;
  char c[16], fpath[256], fname[5], density[16], probe_radius[16];
  char fname_tp[256];
  char c1[10], c2[10], c3[10], c4[10], c5[10];
  double a1, a2, a3, b1, b2, b3;
  double epsp, epsw, bulk_strength, theta;
  int order, maxparnode;
  int mesh_flag,ierr,i;
  double *t_chrpos, *t_atmchr, *t_atmrad;

  extern int tabipb();
/********************************************************/

  fp=fopen("usrdata.in","r");
    ierr=fscanf(fp,"%s %s",c,fname);
    ierr=fscanf(fp,"%s %s",c,density);
    ierr=fscanf(fp,"%s %s",c,probe_radius);
    ierr=fscanf(fp,"%s %lf",c,&epsp);
    ierr=fscanf(fp,"%s %lf",c,&epsw);
    ierr=fscanf(fp,"%s %lf",c,&bulk_strength);
    ierr=fscanf(fp,"%s %d",c,&order);
    ierr=fscanf(fp,"%s %d",c,&maxparnode);
    ierr=fscanf(fp,"%s %lf",c,&theta);
  fclose(fp);

  mesh_flag=1;

/********************************************************/

  sprintf(fpath, "");

  sprintf(fname_tp, "%s%s.pqr", fpath, fname);
  fp = fopen(fname_tp, "r");
  sprintf(fname_tp, "%s%s.xyzr", fpath, fname);
  wfp=fopen(fname_tp, "w");
  int ch, number_of_lines = 0;

  while(fscanf(fp, "%s %s %s %s %s %lf %lf %lf %lf %lf",
               c1, c2, c3, c4, c5, &a1, &a2, &a3, &b1, &b2) != EOF) {
    fprintf(wfp, "%f %f %f %f\n", a1, a2, a3, b2);
    number_of_lines++;
  }

  fclose(fp);
  fclose(wfp);

  sprintf(fname_tp, "%s%s.pqr", fpath, fname);
  fp = fopen(fname_tp, "r");

  if ((t_chrpos = (double *) malloc(3 * number_of_lines * sizeof(double))) == NULL) {
          printf("Error in allocating t_chrpos!\n");
  }
  if ((t_atmchr = (double *) malloc(number_of_lines * sizeof(double))) == NULL) {
          printf("Error in allocating t_atmchr!\n");
  }
  if ((t_atmrad = (double *) malloc(number_of_lines * sizeof(double))) == NULL) {
          printf("Error in allocating t_atmrad!\n");
  }

  for (i = 0; i < number_of_lines; i++) {
          ierr = fscanf(fp, "%s %s %s %s %s %lf %lf %lf %lf %lf",
                        c1,c2,c3,c4,c5,&a1,&a2,&a3,&b1,&b2);
          t_chrpos[3*i] = a1;
          t_chrpos[3*i + 1] = a2;
          t_chrpos[3*i + 2] = a3;
          t_atmchr[i] = b1;
          t_atmrad[i] = b2;
  }

  fclose(fp);
  printf("Finished reading charge (.pqr) file...\n");

  tabipb(fpath, fname, number_of_lines, density, probe_radius, epsp, epsw,
         bulk_strength, order, maxparnode, theta, mesh_flag,
         t_chrpos,t_atmchr,t_atmrad);

  free(t_atmchr);
  free(t_chrpos);
  free(t_atmrad);

  return 0;
}
