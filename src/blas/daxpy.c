/*                                                                                           
 * This file is part of CLAPACK, the f2c'ed version of LAPACK                                
 * hosted at netlib.org.                                                                     
 *                                                                                           
 * For more information, see the LAPACK User's Guide:                                        
 *                                                                                           
  @BOOK{laug,                                                                                
        AUTHOR = {Anderson, E. and Bai, Z. and Bischof, C. and                               
                  Blackford, S. and Demmel, J. and Dongarra, J. and                          
                  Du Croz, J. and Greenbaum, A. and Hammarling, S. and                       
                  McKenney, A. and Sorensen, D.},                                            
        TITLE = {{LAPACK} Users' Guide},                                                     
        EDITION = {Third},                                                                   
        PUBLISHER = {Society for Industrial and Applied Mathematics},                        
        YEAR = {1999},                                                                       
        ADDRESS = {Philadelphia, PA},                                                        
        ISBN = {0-89871-447-8 (paperback)}  }                                                
 *                                                                                           
 */ 

/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
 */

#include "f2c.h"

/* Subroutine */ int daxpy_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx, doublereal *dy, integer *incy)
{


    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     constant times a vector plus a vector.   
       uses unrolled loops for increments equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    if (*n <= 0) {
	return 0;
    }
    if (*da == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	DY(iy) += *da * DX(ix);
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1   


          clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
	DY(i) += *da * DX(i);
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 4) {
	DY(i) += *da * DX(i);
	DY(i + 1) += *da * DX(i + 1);
	DY(i + 2) += *da * DX(i + 2);
	DY(i + 3) += *da * DX(i + 3);
/* L50: */
    }
    return 0;
} /* daxpy_ */

