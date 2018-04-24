/**************************************************************************
* FILE NAME: TABIPBstruct.h                                               *
*                                                                         *
* PURPOSE: typedefs for TABIPB parm and vars structs.                     *
*          TABIPBparm contains parameters specifying run information.     *
*          TABIPBvars contains the atomic information, meshed surfaces,   *
*          and potentials/ their normal derivatives after computation.    *
*          Additionally, contains declaration of Sphinx interface routine *
*                                                                         *
* AUTHORS: Leighton Wilson, University of Michigan, Ann Arbor, MI         *
*          Jiahui Chen, Southern Methodist University, Dallas, TX         *
*                                                                         *
* BASED ON PACKAGE ORIGINALLY WRITTEN IN FORTRAN BY:                      *
*          Weihua Geng, Southern Methodist University, Dallas, TX         *
*          Robery Krasny, University of Michigan, Ann Arbor, MI           *
*                                                                         *
**************************************************************************/

#ifndef H_TABIPB_STRUCT_H
#define H_TABIPB_STRUCT_H

typedef struct sTABIPBparm {

    /* mesh program directory */
    char fpath[256];
    
    /* pqr file location */
    char fname[256];
    
    /* mesh settings */
    double density;
    double probe_radius;

    /* physical parameters */
    double temp;
    double epsp;
    double epsw;
    double bulk_strength;
    
    /* set and used locally */
    double eps;
    double kappa;
    double kappa2;

   /* treecode parameters */
    int order;
    int maxparnode;
    double theta;

   /* mesh program: 0 = msms, 1 = NanoShaper SES, 2 = NanoShaper Skin */
    int mesh_flag;

   /* number of atoms */
    int number_of_lines;

   /* output of potential data */
   /* char 1: DAT, char 2: VTK, char 3: CSV */
    char output_datafile[3];

} TABIPBparm;


typedef struct sTABIPBvars {

    /* solvation and coulombic energy */
    double soleng, couleng;

    /* number of atoms */
    int natm;
    
    /* atomic radii, charges, position */
    double *atmrad, *atmchr, *chrpos;

    /* number of vertices and faces */
    int nspt, nface;

    /* surface area of mesh */
    double surface_area;
    
    /* positions and normals of vertices */
    double **vert, **snrm;

    /* surface potential on vertices */
    double *vert_ptl;
    double max_vert_ptl, min_vert_ptl;
    double max_der_vert_ptl, min_der_vert_ptl;

    /* surface potential on elements */
    double *xvct;
    double max_xvct, min_xvct;
    double max_der_xvct, min_der_xvct;

    /* iterations for GMRes convergence */
    int gmres_iter;

    /* connectivity data for surface triangulation */
    int **face;

} TABIPBvars;


int sphinx2tabipb(TABIPBparm *parm, TABIPBvars *vars);

#endif /* H_TABIPB_STRUCT_H */
