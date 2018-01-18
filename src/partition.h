/**************************************************************************
* FILE NAME: partition.h                                                  *
*                                                                         *
* PURPOSE: header for partition utility used by PartitionEight in         *
*          treecode.c for constructing the tree                           *
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
* 01/15/2018  Leighton Wilson   Created                                   *
*                                                                         *
**************************************************************************/

#ifndef H_PARTITION_H
#define H_PARTITION_H

int Partition(double *a, double *b, double *c, int *indarr,
              int ibeg, int iend, double val);

#endif /* H_PARTITION_H */
