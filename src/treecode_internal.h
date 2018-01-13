/*
 * C header file for treecode internal routines of tabipb
 *
 * C version authored by:
 * Jiahui Chen, Southern Methodist University, Dallas, TX
 * Leighton Wilson, University of Michigan, Ann Arbor, MI
 *
 * Based on package originally written in FORTRAN by:
 * Weihua Geng, Southern Methodist University, Dallas, TX
 * Robery Krasny, University of Michigan, Ann Arbor, MI
 *
 * Last modified by Leighton Wilson, 01/11/2018
 */

#ifndef H_TREECODE_INTERNAL_H
#define H_TREECODE_INTERNAL_H

#include "tree_node_struct.h"


int Setup(double *x, double *y, double *z, double *q, int numpars,
          int order, int iflag, double xyz_limits[6]);

int CreateTree(TreeNode *p, int ibeg, int iend, double xyzmm[6], int level);

int PartitionEight(double xyzmms[6][8], double xl, double yl, double zl,
                   double lmax, double x_mid, double y_mid, double z_mid,
                   int ind[8][2]);

int ComputePBKernel(double *phi);

int ComputeAllMoments(TreeNode *p, int ifirst);

int ComputeMoments(TreeNode *p);

int RunTreecode(TreeNode *p, double peng[2], double *tpoten_old, double tempq[2][16]);

int ComputeTreePB(TreeNode *p, double peng[2], double tempq[2][16]);

int ComputeCoeffs(TreeNode *p, double kappa);

int ComputeDirectPB(double peng[2], int ibeg, int iend, double *tpoten_old);

int RemoveMoments(TreeNode *p);

int RemoveNode(TreeNode *p);


double MinVal(double *variables, int number);

double MaxVal(double *variables, int number);

int Partition(double *a, double *b, double *c, double *q, int *indarr,
              int ibeg, int iend, double val, int numpars);

#endif /* H_TREECODE_INTERNAL_H */
