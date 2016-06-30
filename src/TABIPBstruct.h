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

#ifndef H_TABIPBPARM_H
#define H_TABIPBPARM_H
struct sTABIPBparm {

  /* molecule ID */
  char fpath[256];
  char fname[5];
  char density[16];
  char probe_radius[16];

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

#endif /* H_TABIPBPARM_H */
