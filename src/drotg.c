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

/* Table of constant values */

static doublereal c_b4 = 1.;

/* Subroutine */ int drotg_(doublereal *da, doublereal *db, doublereal *c, 
	doublereal *s)
{


    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal r, scale, z, roe;


/*     construct givens plane rotation.   
       jack dongarra, linpack, 3/11/78. */


    roe = *db;
    if (abs(*da) > abs(*db)) {
	roe = *da;
    }
    scale = abs(*da) + abs(*db);
    if (scale != 0.) {
	goto L10;
    }
    *c = 1.;
    *s = 0.;
    r = 0.;
    z = 0.;
    goto L20;
L10:
/* Computing 2nd power */
    d__1 = *da / scale;
/* Computing 2nd power */
    d__2 = *db / scale;
    r = scale * sqrt(d__1 * d__1 + d__2 * d__2);
    r = d_sign(&c_b4, &roe) * r;
    *c = *da / r;
    *s = *db / r;
    z = 1.;
    if (abs(*da) > abs(*db)) {
	z = *s;
    }
    if (abs(*db) >= abs(*da) && *c != 0.) {
	z = 1. / *c;
    }
L20:
    *da = r;
    *db = z;
    return 0;
} /* drotg_ */

