/*
 * C header file for global variables of tabipb
 *
 * C version authored by:
 * Jiahui Chen, Southern Methodist University, Dallas, TX
 * Leighton Wilson, University of Michigan, Ann Arbor, MI
 *
 * Based on package originally written in FORTRAN by:
 * Weihua Geng, Southern Methodist University, Dallas, TX
 * Robery Krasny, University of Michigan, Ann Arbor, MI
 *
 * Last modified by Leighton Wilson, 06/20/2016
 */

#ifndef H_GLOBAL_PARAMS_H
#define H_GLOBAL_PARAMS_H

/*constant variables */

const double PI = 3.14159265358979324;
const double ONE_OVER_4PI = 0.079577471545948;
const double KCAL_TO_KJ = 4.184;
const double BULK_COEFF = 2529.12179861515279;
const double UNITS_COEFF = 1389.3875744; // 332.0716 * kcal2kj
const double UNITS_PARA = 8729.779593448; // 2 * constUnits * PI

#endif /* H_GLOBAL_PARAMS_H */
