#include <stdio.h>
#include <math.h>

#include "partition.h"

/* 
 * partition determines the index MIDIND, after partitioning in place the arrays a, b, c,
 * and q, such that a(ibeg:midind) <= val and a(midind+1:iend) > val. If on entry, ibeg >
 * iend, or a(ibeg:iend) > val then midind is returned as ibeg-1.
 */


void pc_partition(double *a, double *b, double *c, int *indarr,
                  int ibeg, int iend, double val, int *midind)
{
    double ta, tb, tc, td, tw;
    int lower, upper, tind;
    
    
    if (ibeg < iend) {
        /*
         * temporarily stores ibeg entries and set a(ibeg) = val
         * for the partitioning algorithm.
         */
        ta = a[ibeg];
        tb = b[ibeg];
        tc = c[ibeg];
        tind = indarr[ibeg];
        a[ibeg] = val;
        upper = ibeg;
        lower = iend;
        
        while (upper != lower) {
            while ((upper < lower) && (val < a[lower])) {
                lower--;
            }
            
            if (upper != lower) {
                a[upper] = a[lower];
                b[upper] = b[lower];
                c[upper] = c[lower];
                indarr[upper] = indarr[lower];
            }
            
            while ((upper < lower) && (val >= a[upper])) {
                upper++;
            }
            
            if (upper != lower) {
                a[lower] = a[upper];
                b[lower] = b[upper];
                c[lower] = c[upper];
                indarr[lower] = indarr[upper];
            }
        }
        
        *midind = upper;
        
        
        /* replace TA in position upper and change midind if ta > val */
        
        if (ta > val)
            *midind = upper - 1;
        
        
        a[upper] = ta;
        b[upper] = tb;
        c[upper] = tc;
        indarr[upper] = tind;
        
    } else if (ibeg == iend) {
        
        if (a[ibeg] < val)
            *midind = ibeg;
        else
            *midind = ibeg - 1;
        
    } else {
        
        *midind = ibeg - 1;
        
    }
    
    return;
    
} /* END function pc_partition */


void pc_partition_8(double *x, double *y, double *z, int *orderarr, double xyzmms[6][8],
                    double xl, double yl, double zl, int *numposchild, int max_num_children,
                    double x_mid, double y_mid, double z_mid, int ind[8][2])
{
    int temp_ind;

    *numposchild = 1;
    
    double lmax = xl;
    if (lmax < yl) lmax = yl;
    if (lmax < zl) lmax = zl;

    double critlen = lmax / sqrt(2.0);
    
    int divide_x = 0;
    int divide_y = 0;
    int divide_z = 0;
    
    if (xl >= critlen) divide_x = 1;
    if (yl >= critlen) divide_y = 1;
    if (zl >= critlen) divide_z = 1;
    
    if (max_num_children == 4) {
        if (xl < yl && xl < zl) divide_x = 0;
        if (yl < xl && yl < zl) divide_y = 0;
        if (zl < xl && zl < yl) divide_z = 0;
        
        if (divide_x + divide_y + divide_z == 3) divide_x = 0;
    }
    
    if (max_num_children == 2) {
        if (xl < yl || xl < zl) divide_x = 0;
        if (yl < xl || yl < zl) divide_y = 0;
        if (zl < xl || zl < yl) divide_z = 0;
        
        if (divide_x + divide_y + divide_z == 3) {
            divide_x = 0;
            divide_y = 0;
        }
        
        if (divide_x + divide_y + divide_z == 2) {
            if        (divide_x == 0) {
                divide_y = 0;
            } else if (divide_y == 0) {
                divide_z = 0;
            } else if (divide_z == 0) {
                divide_x = 0;
            }
        }
    }

    if (xl >= critlen) {

        pc_partition(x, y, z, orderarr, ind[0][0], ind[0][1],
                     x_mid, &temp_ind);

        ind[1][0] = temp_ind + 1;
        ind[1][1] = ind[0][1];
        ind[0][1] = temp_ind;

        for (int i = 0; i < 6; i++)
            xyzmms[i][1] = xyzmms[i][0];

        xyzmms[1][0] = x_mid;
        xyzmms[0][1] = x_mid;
        *numposchild = 2 * *numposchild;

    }

    if (yl >= critlen) {

        for (int i = 0; i < *numposchild; i++) {
            pc_partition(y, x, z, orderarr, ind[i][0], ind[i][1],
                         y_mid, &temp_ind);

            ind[*numposchild + i][0] = temp_ind + 1;
            ind[*numposchild + i][1] = ind[i][1];
            ind[i][1] = temp_ind;

            for (int j = 0; j < 6; j++)
                xyzmms[j][*numposchild + i] = xyzmms[j][i];

            xyzmms[3][i] = y_mid;
            xyzmms[2][*numposchild + i] = y_mid;
        }

        *numposchild = 2 * *numposchild;

    }

    if (zl >= critlen) {

        for (int i = 0; i < *numposchild; i++) {
            pc_partition(z, x, y, orderarr, ind[i][0], ind[i][1],
                         z_mid, &temp_ind);

            ind[*numposchild + i][0] = temp_ind + 1;
            ind[*numposchild + i][1] = ind[i][1];
            ind[i][1] = temp_ind;

            for (int j = 0; j < 6; j++)
                xyzmms[j][*numposchild + i] = xyzmms[j][i];

            xyzmms[5][i] = z_mid;
            xyzmms[4][*numposchild + i] = z_mid;
        }

        *numposchild = 2 * *numposchild;

    }

    return;

} /* END of function partition_8 */
