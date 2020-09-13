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


#include "f2c.h"

#ifdef KR_headers
double d_sign(a,b) doublereal *a, *b;
#else
double d_sign(doublereal *a, doublereal *b)
#endif
{
double x;
x = (*a >= 0 ? *a : - *a);
return( *b >= 0 ? x : -x);
}
