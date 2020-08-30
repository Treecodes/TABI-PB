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

#include <vector>

struct Particles {
    
    int num;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    
    std::vector<double> nx;
    std::vector<double> ny;
    std::vector<double> nz;
    
    std::vector<double> area;
    std::vector<double> source_term;
    std::vector<double> xvct;

    std::vector<size_t> order;
};

#endif /* H_PARTICLE_STRUCT_H */
