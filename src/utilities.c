/*
 * C helper routines for treecode in tabipb
 *
 * C version authored by:
 * Leighton Wilson, University of Michigan, Ann Arbor, MI
 * Jiahui Chen, Southern Methodist University, Dallas, TX
 *
 * Based on package originally written in FORTRAN by:
 * Weihua Geng, Southern Methodist University, Dallas, TX
 * Robery Krasny, University of Michigan, Ann Arbor, MI
 *
 * Added TriangleArea by Leighton Wilson, 01/15/2018
 * Created by Leighton Wilson, 01/12/2018
 */

#include <math.h>

#include "utilities.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* helper functions                                          * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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

    return(area);
}
/**********************************************************/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
