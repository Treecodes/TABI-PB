/**************************************************************************
* FILE NAME: xmalloc.c                                                    *
*                                                                         *
* PURPOSE: safe memory allocation routine used by array.h                 *
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

#include <stdio.h>
#include "xmalloc.h"

void *malloc_or_exit(size_t nbytes, const char *file, int line)
{
    void *x;
    if ((x = malloc(nbytes)) == NULL || nbytes == 0) {
        fprintf(stderr, "%s: line %d: malloc() of %zu bytes failed\n",
                file, line, nbytes);
        exit(EXIT_FAILURE);
    } else {
        return x;
    }
}
