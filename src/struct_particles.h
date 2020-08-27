/**************************************************************************
* FILE NAME: particle_struct.h                                            *
*                                                                         *
* PURPOSE: header file for tree particle struct, used to store particle   *
*          information for treecode and interface with GMRes              *
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
* 01/12/2018  Leighton Wilson   Renamed struct, added source term         *
* 01/10/2018  Leighton Wilson   Created                                   *
*                                                                         *
**************************************************************************/

#ifndef H_PARTICLES_STRUCT_H
#define H_PARTICLES_STRUCT_H

struct Particles {
    
    int num;
    double *x;
    double *y;
    double *z;
    
    double *nx;
    double *ny;
    double *nz;
    
    double *area;
    double *source_term;
    double *xvct;

    int *order;
};

#endif /* H_PARTICLE_STRUCT_H */
