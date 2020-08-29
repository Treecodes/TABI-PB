/**************************************************************************
* FILE NAME: utilities.c                                                  *
*                                                                         *
* PURPOSE: Contains helper routines for calculating min and max values    *
*          in an array (used by treecode.c and tabipb.c), as well as      *
*          routine to calculate triangle area (used by readin.c and       *
*          tabipb.c)                                                      *
*                                                                         *
* AUTHORS: Leighton Wilson, University of Michigan, Ann Arbor, MI         *
*          Jiahui Chen, Southern Methodist University, Dallas, TX         *
*                                                                         *
* BASED ON PACKAGE ORIGINALLY WRITTEN IN FORTRAN BY:                      *
*          Weihua Geng, Southern Methodist University, Dallas, TX         *
*          Robery Krasny, University of Michigan, Ann Arbor, MI           *
*                                                                         *
**************************************************************************/

#include <math.h>

#include "utilities.h"

/**********************************************************/
double MinVal(double *variables, int number) 
{
    int i;
    double min_val;

    min_val = variables[0];
    for (i = 1; i < number; i++) {
        if (min_val > variables[i]) {
            min_val = variables[i];
        }
    }
    
    return min_val;
}
/**********************************************************/


/**********************************************************/
double MaxVal(double *variables, int number) 
{
    int i;
    double max_val;
  
    max_val = variables[0];
    for (i = 1; i < number; i++) {
        if (max_val < variables[i])
            max_val = variables[i];
    }
    
    return max_val;
}
/**********************************************************/


/**********************************************************/
/* function computing the area of a triangle given vertices coodinates */
double TriangleArea(double v[3][3])
{
    int i;
    double a[3], b[3], c[3], aa, bb, cc, ss, area;

    for (i = 0; i <= 2; i++) {
        a[i] = v[i][0] - v[i][1];
        b[i] = v[i][0] - v[i][2];
        c[i] = v[i][1] - v[i][2];
    }

    aa = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    bb = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
    cc = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);

    ss = 0.5 * (aa + bb + cc);
    area = sqrt(ss * (ss-aa) * (ss-bb) * (ss-cc));

    return area;
}
/**********************************************************/
