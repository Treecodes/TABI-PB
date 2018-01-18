/**************************************************************************
* FILE NAME: tabipb.h                                                     *
*                                                                         *
* PURPOSE: header for calling TABI-PB solver                              *
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
* 01/11/2018  Leighton Wilson   Created                                   *
*                                                                         *
**************************************************************************/

#ifndef H_TABIPB_H
#define H_TABIPB_H

#include "TABIPBstruct.h"

int TABIPB(TABIPBparm *parm, TABIPBvars *vars);

#endif /* H_TABIPB_H */
