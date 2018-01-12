/*
 * C version authored by:
 * Jiahui Chen, Southern Methodist University, Dallas, TX
 * Leighton Wilson, University of Michigan, Ann Arbor, MI
 *
 * Based on package originally written in FORTRAN by:
 * Weihua Geng, Southern Methodist University, Dallas, TX
 * Robery Krasny, University of Michigan, Ann Arbor, MI
 *
 * Works for Sphinx by Jiahui at 7/14/2016
 * Rebuild the architecture of wrapper by Jiahui at 6/30/2016
 * Build matrix free and nanoshaper by Leighton at 6/23/2016
 *
 */

#include <string.h>
#include "TABIPBstruct.h"

int main(int argc, char **argv)
{
  /* main reads the input file, writes xyzr file for msms and sets up position,
     raduis and charges */

    FILE *fp, *wfp;
    char c[16];
    char fname_tp[256];
    char c1[10], c2[10], c3[10], c4[10], c5[10];
    double a1, a2, a3, b1, b2, b3;
    double density, radius, epsp, epsw, bulk_strength, theta, temp;
    int maxparnode, order, mesh_flag, output_datafile, ierr, i;

  /* timing functions for *nix systems */
#ifndef _WIN32                                                                     
    extern void timer_start();
    extern void timer_end();
#endif  

    extern int OutputPrint();
    extern int OutputVTK();

#ifndef _WIN32                                                                     
    timer_start("TOTAL_TIME");
#endif  

    TABIPBparm *main_parm = malloc(sizeof *main_parm);
    TABIPBvars *main_vars = malloc(sizeof *main_vars);

    extern int tabipb();
/********************************************************/

    fp = fopen("usrdata.in", "r");
    ierr = fscanf(fp, "%s %s", c, main_parm->fname);

    ierr = fscanf(fp, "%s %lf", c, &density);
    main_parm->density = density;

    ierr = fscanf(fp, "%s %lf", c, &radius);
    main_parm->probe_radius = radius;

    ierr = fscanf(fp, "%s %lf", c, &epsp);
    main_parm->epsp = epsp;

    ierr = fscanf(fp, "%s %lf", c, &epsw);
    main_parm->epsw = epsw;

    ierr = fscanf(fp, "%s %lf", c, &bulk_strength);
    main_parm->bulk_strength = bulk_strength;

    ierr = fscanf(fp, "%s %lf", c, &temp);
    main_parm->temp = temp;

    ierr = fscanf(fp, "%s %d", c, &order);
    main_parm->order = order;

    ierr = fscanf(fp, "%s %d", c, &maxparnode);
    main_parm->maxparnode = maxparnode;

    ierr = fscanf(fp, "%s %lf", c, &theta);
    main_parm->theta = theta;

    ierr = fscanf(fp, "%s %d", c, &mesh_flag);
    main_parm->mesh_flag = mesh_flag;

    ierr = fscanf(fp, "%s %d", c, &output_datafile);
    main_parm->output_datafile = output_datafile;
    fclose(fp);

/********************************************************/
    sprintf(main_parm->fpath, "");

    sprintf(fname_tp, "%s%s.pqr", main_parm->fpath, main_parm->fname);
    fp = fopen(fname_tp, "r");

    sprintf(fname_tp, "molecule.xyzr");
    wfp = fopen(fname_tp, "w");
    
    int ch;// main_parm->number_of_lines = 0;

    while (fscanf(fp, "%s %s %s %s %s %lf %lf %lf %lf %lf",
           c1, c2, c3, c4, c5, &a1, &a2, &a3, &b1, &b2) != EOF) {
        if (strncmp(c1, "ATOM", 4) == 0) {
            fprintf(wfp, "%f %f %f %f\n", a1, a2, a3, b2);
            main_parm->number_of_lines++;
        }
    }

    fclose(fp);
    fclose(wfp);
    printf("Finished assembling atomic information (.xyzr) file...\n");

    sprintf(fname_tp, "%s%s.pqr", main_parm->fpath, main_parm->fname);
    fp = fopen(fname_tp, "r");

    if ((main_vars->chrpos = (double *) calloc(3 * main_parm->number_of_lines, sizeof(double))) == NULL) {
        printf("Error in allocating t_chrpos!\n");
    }
    if ((main_vars->atmchr = (double *) calloc(main_parm->number_of_lines, sizeof(double))) == NULL) {
        printf("Error in allocating t_atmchr!\n");
    }
    if ((main_vars->atmrad = (double *) calloc(main_parm->number_of_lines, sizeof(double))) == NULL) {
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
    printf("Finished assembling charge structures from .pqr file...\n");

    ierr = tabipb(main_parm, main_vars);

    ierr = OutputPrint(main_vars);
    if (output_datafile == 1) {
        ierr=OutputVTK(main_parm, main_vars);
    }

    free(main_parm);
    free(main_vars->atmchr);
    free(main_vars->chrpos);
    free(main_vars->atmrad);
    free(main_vars->vert_ptl); // allocate in output_potential()
    free(main_vars->xvct);
    free_matrix(main_vars->vert);
    free_matrix(main_vars->snrm);
    free_matrix(main_vars->face);
    free(main_vars);

#ifndef _WIN32
    timer_end();
#endif

    return 0;
}
