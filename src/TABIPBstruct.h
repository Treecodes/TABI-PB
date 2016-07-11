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
 * Rebuild the architecture of wrapper at 6/30/2016
 */

#ifndef H_TABIPBSTRUCT_H
#define H_TABIPBSTRUCT_H

#include <stdlib.h> /* calloc() */
#include <stdio.h> /* FILE */
#include "array.h"

struct sTABIPBparm {

  /* molecule ID */
  char fpath[256];
  char fname[5];
  char density[16];
  char probe_radius[16];

  double temp;
  double epsp;
  double epsw;
  double bulk_strength;

  /* treecode_parm */
  int order;
  int maxparnode;
  double theta;

  /* 0 msms and 1 NanoShaper */
  int mesh_flag;

  int number_of_lines;

};

typedef struct sTABIPBparm TABIPBparm;

struct sTABIPBvars {

  /* solvation energy */
  double soleng;

  /* number of nodes, number of triangles */
  int nspt, nface, natm;

  /* msms variables */
  int **extr_v, **extr_f; //[3][nspt],[2][nface]

  /*
  double *tr_xyz, *tr_q; //[3]*[nface]
  double *tr_area, *bvct;
   */
  double *atmrad, *atmchr, *chrpos;

  /* vertices, normal vector*/
  double **vert, **snrm;

  /* surface potential on vertices area */
  double *vert_ptl;
  double max_vert_ptl, min_vert_ptl;
  double max_der_vert_ptl, min_der_vert_ptl;

  /* surface potential on elements area */
  double *xvct;
  double max_xvct, min_xvct;
  double max_der_xvct, min_der_xvct;

  /* connectivity data for MSMS surface triangulation */
  int **face;

};

typedef struct sTABIPBvars TABIPBvars;

#endif /* H_TABIPBSTRUCT_H */
