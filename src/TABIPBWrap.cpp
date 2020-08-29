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
**************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "readin.h"
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

    sprintf(fname_tp, "molecule.xyzr");
    wfp = fopen(fname_tp,"w");
    for (i = 0; i < parm->number_of_lines; i++) {
        fprintf(wfp, "%f %f %f %f\n", vars->chrpos[3*i], vars->chrpos[3*i + 1],
                vars->chrpos[3*i + 2], vars->atmrad[i]);
    }
    fclose(wfp);
    
    ierr = Readin(parm, vars);
    
    ierr = TABIPB(parm, vars);

    ierr = OutputPrint(vars);
    
    /* REMEMBER: These comparisons are now characters */
    
    if (parm->output_datafile[0] == '1') {
        ierr = OutputDAT("tabipb_output", vars);
    } else if (parm->output_datafile[1] == '1') {
        ierr = OutputVTK("tabipb_output", vars);
    }

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

    free(vars->vert[0]);
    free(vars->vert[1]);
    free(vars->vert[2]);

    free(vars->snrm[0]);
    free(vars->snrm[1]);
    free(vars->snrm[2]);

    free(vars->face[0]);
    free(vars->face[1]);
    free(vars->face[2]);

    return 0;
}
/**********************************************************/
