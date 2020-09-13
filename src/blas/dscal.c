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

/* Subroutine */ int dscal_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx)
{


    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, m, nincx, mp1;


/*     scales a vector by a constant.   
       uses unrolled loops for increment equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 3/93 to return if incx .le. 0.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DX(I) dx[(I)-1]


    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {
	DX(i) = *da * DX(i);
/* L10: */
    }
    return 0;

/*        code for increment equal to 1   


          clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i = 1; i <= m; ++i) {
	DX(i) = *da * DX(i);
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i = mp1; i <= *n; i += 5) {
	DX(i) = *da * DX(i);
	DX(i + 1) = *da * DX(i + 1);
	DX(i + 2) = *da * DX(i + 2);
	DX(i + 3) = *da * DX(i + 3);
	DX(i + 4) = *da * DX(i + 4);
/* L50: */
    }
    return 0;
} /* dscal_ */

