/**************************************************************************
* FILE NAME: global_params.h                                              *
*                                                                         *
* PURPOSE: global parameters used to define units, pi                     *
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
* 01/12/2018  Leighton Wilson   Created, moved from tabipb header         *
*                                                                         *
**************************************************************************/

#ifndef H_GLOBAL_PARAMS_H
#define H_GLOBAL_PARAMS_H

#define PI 3.14159265358979324
#define ONE_OVER_4PI 0.079577471545948
#define KCAL_TO_KJ 4.184
#define BULK_COEFF 2529.12179861515279
#define UNITS_COEFF 1389.3875744 // 332.0716 * kcal2kj
#define UNITS_PARA 8729.779593448 // 2 * constUnits * PI

#endif /* H_GLOBAL_PARAMS_H */
