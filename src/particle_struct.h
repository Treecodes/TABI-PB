/*
 * C header file for tree particle struct
 *
 * C version authored by:
 * Jiahui Chen, Southern Methodist University, Dallas, TX
 * Leighton Wilson, University of Michigan, Ann Arbor, MI
 *
 * Based on package originally written in FORTRAN by:
 * Weihua Geng, Southern Methodist University, Dallas, TX
 * Robery Krasny, University of Michigan, Ann Arbor, MI
 *
 * Last modified by Leighton Wilson, 01/11/2018
 * Created by Leighton Wilson, 01/10/2018
 */

#ifndef H_PARTICLE_STRUCT_H
#define H_PARTICLE_STRUCT_H


typedef struct sTreeParticles {
    
    int num_particles;
    double **particle_position;
    double **particle_normal;
    double *particle_area;
    double *source_term;

} TreeParticles;

#endif /* H_PARTICLE_STRUCT_H */
