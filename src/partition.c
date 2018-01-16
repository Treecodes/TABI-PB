/*
 * C partition routine for treecode in tabipb
 *
 * C version authored by:
 * Leighton Wilson, University of Michigan, Ann Arbor, MI
 * Jiahui Chen, Southern Methodist University, Dallas, TX
 *
 * Based on package originally written in FORTRAN by:
 * Weihua Geng, Southern Methodist University, Dallas, TX
 * Robery Krasny, University of Michigan, Ann Arbor, MI
 *
 * Created by Leighton Wilson, 01/15/2018
 */

#include "partition.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* partition function                                        * * * */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/********************************************************/
int Partition(double *a, double *b, double *c, int *indarr,
              int ibeg, int iend, double val)
{
/* PARTITION determines the index MIDIND, after partitioning
 * in place the  arrays A,B,C and Q,  such that
 * A(IBEG:MIDIND) <= VAL and  A(MIDIND+1:IEND) > VAL.
 * If on entry IBEG > IEND or  A(IBEG:IEND) > VAL then MIDIND
 * is returned as IBEG-1.  */
    double ta, tb, tc;
    int lower, upper, tind;
    int midind;

    if (ibeg < iend) {
/* temporarily store IBEG entries and set A(IBEG)=VAL for
 * the partitoning algorithm.  */
        ta = a[ibeg];
        tb = b[ibeg];
        tc = c[ibeg];
        tind = indarr[ibeg];
        a[ibeg] = val;/*val=mid val on that direction*/
        upper = ibeg;
        lower = iend;

        while (upper != lower) {
            while (upper < lower && val < a[lower])
                lower -= 1;
            if (upper != lower) {
                a[upper] = a[lower];
                b[upper] = b[lower];
                c[upper] = c[lower];
                indarr[upper] = indarr[lower];
            }
            while (upper < lower && val >= a[upper])
            upper += 1;
            if (upper != lower) {
                a[lower] = a[upper];
                b[lower] = b[upper];
                c[lower] = c[upper];
                indarr[lower] = indarr[upper];
            }
        }
        midind = upper;
/* replace TA in position UPPER and change MIDIND if TA > VAL */
        if (ta > val) {
            midind = upper - 1;
        }

        a[upper] = ta;
        b[upper] = tb;
        c[upper] = tc;
        indarr[upper] = tind;
    
    } else if (ibeg == iend) {
    
        if (a[ibeg] <= val) {
            midind = ibeg;
        } else {
            midind = ibeg - 1;
        }
    
    } else {
        midind = ibeg - 1;
    }

    return (midind);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
