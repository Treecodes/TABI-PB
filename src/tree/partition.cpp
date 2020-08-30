#include <stdio.h>
#include <math.h>
#include <iostream>
#include <numeric>

#include "partition.h"

/* 
 * partition determines the index MIDIND, after partitioning in place the arrays a, b, c,
 * and q, such that a(ibeg:midind) <= val and a(midind+1:iend) > val. If on entry, ibeg >
 * iend, or a(ibeg:iend) > val then midind is returned as ibeg-1.
 */
 
 
template< typename order_iterator, typename value_iterator >
void reorder_inplace_destructive( order_iterator order_begin, order_iterator order_end,
     order_iterator i, value_iterator v, value_iterator w, value_iterator x)// value_iterator w, value_iterator x)
{
    typedef typename std::iterator_traits< value_iterator >::value_type value_t;
    typedef typename std::iterator_traits< order_iterator >::value_type index_t;
    typedef typename std::iterator_traits< order_iterator >::difference_type diff_t;
    
    diff_t remaining = order_end - 1 - order_begin;
    for ( index_t s = index_t(); remaining > 0; ++ s ) {
        index_t d = order_begin[s];
        if ( d == (diff_t) -1 ) continue;
        -- remaining;
        value_t tempv = v[s];
        value_t tempw = w[s];
        value_t tempx = x[s];
        index_t tempi = i[s];
        for ( index_t d2; d != s; d = d2 ) {
            std::swap( tempv, v[d] );
            std::swap( tempw, w[d] );
            std::swap( tempx, x[d] );
            std::swap( tempi, i[d] );
            std::swap( order_begin[d], d2 = (diff_t) -1 );
            -- remaining;
        }
        v[s] = tempv;
        w[s] = tempw;
        x[s] = tempx;
        i[s] = tempi;
    }
}


void pc_partition(double *a, double *b, double *c, size_t *indarr,
                  int ibeg, int iend, double val, int *midind)
{
    double ta, tb, tc;
    size_t lower, upper, tind;
    
    
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


void pc_partition_8(std::vector<double> &x, std::vector <double> &y, std::vector<double> &z,
                    std::vector<size_t> &orderarr, double xyzmms[6][8],
                    double xl, double yl, double zl, int *numposchild,
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

    if (divide_x) {

//  This, unfortunately, does not quite work, but it should be something like this to reorder them in an STL way
//
//        std::vector<size_t> reorder_vec(ind[0][1] - ind[0][0] + 1);
//        std::iota(reorder_vec.begin(), reorder_vec.end(), ind[0][0]);
//
//        auto pivot = std::partition(reorder_vec.begin(), reorder_vec.end(),
//                                    [&x, &x_mid, &ind](size_t elem){ return x[elem] < x_mid; });
//
//        std::transform(reorder_vec.begin(),reorder_vec.end(),reorder_vec.begin(),
//                       [&ind](size_t i){ return i-ind[0][0]; });
//
//        reorder_inplace_destructive(reorder_vec.begin(), reorder_vec.end(),
//                        orderarr.begin()+ind[0][0], x.begin()+ind[0][0], y.begin()+ind[0][0], z.begin()+ind[0][0]);
//
//        ind[1][0] = *pivot + ind[0][0];
//        ind[1][1] = ind[0][1];
//        ind[0][1] = *pivot + ind[0][0] - 1;


        pc_partition(x.data(), y.data(), z.data(), orderarr.data(), ind[0][0], ind[0][1],
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

    if (divide_y) {

        for (int i = 0; i < *numposchild; i++) {
        
            pc_partition(y.data(), x.data(), z.data(), orderarr.data(), ind[i][0], ind[i][1],
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

    if (divide_z) {

        for (int i = 0; i < *numposchild; i++) {
        
            pc_partition(z.data(), x.data(), y.data(), orderarr.data(), ind[i][0], ind[i][1],
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
