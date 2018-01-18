/**************************************************************************
* FILE NAME: TABIPBWrap.c                                                 *
*                                                                         *
* PURPOSE: routines to interface TABI-PB caller with APBS and Sphinx      *
*                                                                         *
* AUTHORS: Leighton Wilson, University of Michigan, Ann Arbor, MI         *
*          Jiahui Chen, Southern Methodist University, Dallas, TX         *
*                                                                         *
* BASED ON PACKAGE ORIGINALLY WRITTEN IN FORTRAN BY:                      *
*          Weihua Geng, Southern Methodist University, Dallas, TX         *
*          Robery Krasny, University of Michigan, Ann Arbor, MI           *
*                                                                         *
* DEVELOPMENT HISTORY:                                                    *
*                                                                         *
* Date        Author            Description Of Change                     *
* ----        ------            ---------------------                     *
* 01/12/2018  Leighton Wilson   Modified included headers                 *
* 06/29/2016  Jiahui Chen       Rebuilding wrapper architecture           *
*                                                                         *
**************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "tabipb.h"
#include "print_output.h"

#include "array.h"
#include "TABIPBstruct.h"

/**********************************************************/
int apbs2tabipb_(TABIPBparm *parm, TABIPBvars *vars)
{
    FILE *wfp;
    char fname_tp[256];
    int i, ierr;

    //sprintf(fname_tp, "%s%s.xyzr",parm->fpath, parm->fname);
    sprintf(fname_tp, "molecule.xyzr");
    wfp = fopen(fname_tp,"w");
    for (i = 0; i < parm->number_of_lines; i++) {
        fprintf(wfp, "%f %f %f %f\n", vars->chrpos[3*i], vars->chrpos[3*i + 1],
                vars->chrpos[3*i + 2], vars->atmrad[i]);
    }
    fclose(wfp);

    ierr = TABIPB(parm, vars);

    ierr = OutputPrint(vars);
    if (parm->output_datafile == 1) ierr = OutputVTK(parm, vars);

    return 0;
}
/**********************************************************/


/**********************************************************/
int sphinx2tabipb(TABIPBparm *parm, TABIPBvars *vars)
{
    FILE *wfp;
    char fname_tp[256];
    int i, ierr;

    sprintf(fname_tp, "molecule.xyzr");
    wfp = fopen(fname_tp,"w");
    for (i = 0; i < parm->number_of_lines; i++) {
        fprintf(wfp, "%f %f %f %f\n", vars->chrpos[3*i], vars->chrpos[3*i + 1],
                vars->chrpos[3*i + 2], vars->atmrad[i]);
    }
    fclose(wfp);

    ierr = TABIPB(parm, vars);

    free(vars->vert_ptl);
    free(vars->xvct);
    free_matrix(vars->vert);
    free_matrix(vars->snrm);
    free_matrix(vars->face);

    return 0;
}
/**********************************************************/
