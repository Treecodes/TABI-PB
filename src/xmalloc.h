/**************************************************************************
* FILE NAME: xmalloc.h                                                    *
*                                                                         *
* PURPOSE: header for xmalloc.c memory allocation routine                 *
*                                                                         *
* AUTHORS: Leighton Wilson, University of Michigan, Ann Arbor, MI         *
*          Jiahui Chen, Southern Methodist University, Dallas, TX         *
*                                                                         *
* BASED ON WORK OF:                                                       *
*          Rouben Rostamian, presented in "Programming Projects in C for  *
*          Students of Engineering, Science, and Mathematics"             *
*                                                                         *
* DEVELOPMENT HISTORY:                                                    *
*                                                                         *
* Date        Author            Description Of Change                     *
* ----        ------            ---------------------                     *
* 06/23/2016  Leighton Wilson   Created                                   *
*                                                                         *
**************************************************************************/

#ifndef H_XMALLOC_H
#define H_XMALLOC_H

#include <stdlib.h>

void *malloc_or_exit(size_t nbytes, const char *file, int line);
#define xmalloc(nbytes)  malloc_or_exit((nbytes), __FILE__, __LINE__)

#endif /*H_XMALLOC_H*/
