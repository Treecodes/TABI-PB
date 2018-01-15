/*
 * C routines for printing tabipb output
 *
 * C version authored by:
 * Leighton Wilson, University of Michigan, Ann Arbor, MI
 * Jiahui Chen, Southern Methodist University, Dallas, TX
 *
 * Based on package originally written in FORTRAN by:
 * Weihua Geng, Southern Methodist University, Dallas, TX
 * Robery Krasny, University of Michigan, Ann Arbor, MI
 *
 * Last modified by Leighton Wilson, 01/12/2018
 */

#include <stdio.h>

#include "print_output.h"

#include "global_params.h"
#include "TABIPBstruct.h"


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* output printing functions                                 * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/************************************/
int OutputPrint(TABIPBvars *vars)
{
    int i;

    printf("\nSolvation energy = %f kJ/mol", vars->soleng);
    printf("\nFree energy = %f kJ/mol\n\n", vars->soleng+vars->couleng);
    printf("The max and min potential and normal derivatives on elements area:\n");
    printf("potential %f %f\n", vars->max_xvct, vars->min_xvct);
    printf("norm derv %f %f\n\n", vars->max_der_xvct,
                                  vars->min_der_xvct);
    printf("The max and min potential and normal derivatives on vertices area:\n");
    printf("potential %f %f\n", vars->max_vert_ptl, vars->min_vert_ptl);
    printf("norm derv %f %f\n\n", vars->max_der_vert_ptl,
                                  vars->min_der_vert_ptl);

    FILE *fp = fopen("surface_potential.dat", "w");
    fprintf(fp, "%d %d\n", vars->nspt, vars->nface);

    for (i = 0; i < vars->nspt; i++)
        fprintf(fp, "%d %f %f %f %f %f %f %f %f\n", i,
                vars->vert[0][i], vars->vert[1][i], vars->vert[2][i],
                vars->snrm[0][i], vars->snrm[1][i], vars->snrm[2][i],
                vars->vert_ptl[i], vars->vert_ptl[i + vars->nspt]);

    for (i = 0; i < vars->nface; i++)
        fprintf(fp, "%d %d %d\n", vars->face[0][i], vars->face[1][i], vars->face[2][i]);

    fclose(fp);

    return 0;
}


/************************************/
int OutputVTK(TABIPBparm *parm, TABIPBvars *vars)
{
    char i_char1[20], nspt_str[20],
         nface_str[20], nface4_str[20];
    int i;

    sprintf(nspt_str, "%d", vars->nspt);
    sprintf(nface_str, "%d", vars->nface);
    sprintf(nface4_str, "%d", vars->nface * 4);

    sprintf(i_char1, "mesh flag: %d", parm->mesh_flag);

    FILE *fp = fopen("surface_potential.vtk", "w");

    fprintf(fp, "# vtk DataFile Version 1.0\n");
    fprintf(fp, "mesh for protein %s, with %s\n", parm->fname, i_char1);
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET POLYDATA\n\n");

    fprintf(fp, "POINTS %s double\n", nspt_str);
    for (i = 0; i < vars->nspt; i++) {
        fprintf(fp, "%f %f %f\n", vars->vert[0][i], vars->vert[1][i],
                                  vars->vert[2][i]);
    }

    fprintf(fp, "POLYGONS %s %s\n", nface_str, nface4_str);
    for (i = 0; i < vars->nface; i++) {
        fprintf(fp, "3 %d %d %d\n", vars->face[0][i] - 1, vars->face[1][i] - 1,
                                    vars->face[2][i] - 1);
    }

    fprintf(fp, "\nPOINT_DATA %s\n", nspt_str);
    fprintf(fp, "SCALARS PotentialVert double\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (i = 0; i < vars->nspt; i++) {
         fprintf(fp, "%f\n", KCAL_TO_KJ * vars->vert_ptl[i]);
    }

    fprintf(fp, "SCALARS NormalPotentialVert double\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (i = 0; i < vars->nspt; i++) {
        fprintf(fp, "%f\n", KCAL_TO_KJ * vars->vert_ptl[vars->nspt + i]);
    }

    //if we want induced surface charges, we can multiply vertnorm by (1/eps + 1)
    fprintf(fp, "\nNORMALS VertNorms double\n");
    for (i = 0; i < vars->nspt; i++) {
        fprintf(fp, "%f %f %f\n", vars->snrm[0][i], vars->snrm[1][i],
                                  vars->snrm[2][i]);
    }

    fprintf(fp, "\nCELL_DATA %s\n", nface_str);
    fprintf(fp, "SCALARS PotentialFace double\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (i = 0; i < vars->nface; i++) {
        fprintf(fp, "%f\n", KCAL_TO_KJ * vars->xvct[i]);
    }

    //if we want induced surface charges, we can multiply vertnorm by (1/eps + 1)
    fprintf(fp, "SCALARS NormalPotentialFace double\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (i = 0; i < vars->nface; i++) {
        fprintf(fp, "%f\n", KCAL_TO_KJ * vars->xvct[vars->nface + i]);
    }

    fclose(fp);

    return 0;
}