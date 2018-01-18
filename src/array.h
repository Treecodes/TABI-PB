/**************************************************************************
* FILE NAME: array.h                                                      *
*                                                                         *
* PURPOSE: header containing macros for vector and array creation using   *
*          the xmalloc.c memory allocation routine                        *
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

#ifndef H_ARRAY_H
#define H_ARRAY_H
#include "xmalloc.h"

#define make_vector(v,n) ((v) = xmalloc((n) * sizeof *(v)))

#define free_vector(v)  do { free(v); v = NULL; } while (0)

#define make_matrix(a, m, n) do {                              \
    size_t make_matrix_loop_counter;                           \
    make_vector(a, (m) + 1);                                   \
    for (make_matrix_loop_counter = 0;                         \
            make_matrix_loop_counter < (m);                    \
            make_matrix_loop_counter++)                        \
        make_vector((a)[make_matrix_loop_counter], (n));       \
    (a)[m] = NULL;                                             \
} while (0)

#define free_matrix(a) do {                                    \
    if (a != NULL) {                                           \
        size_t make_matrix_loop_counter;                       \
        for (make_matrix_loop_counter = 0;                     \
                (a)[make_matrix_loop_counter] != NULL;         \
                make_matrix_loop_counter++)                    \
            free_vector((a)[make_matrix_loop_counter]);        \
        free_vector(a);                                        \
        a = NULL;                                              \
    }                                                          \
} while (0)

#define make_3array(a, l, m, n) do {                           \
    size_t make_3array_loop_counter;                           \
    make_vector(a, (l) + 1);                                   \
    for (make_3array_loop_counter = 0;                         \
            make_3array_loop_counter < (l);                    \
            make_3array_loop_counter++)                        \
        make_matrix((a)[make_3array_loop_counter], (m), (n));  \
    (a)[l] = NULL;                                             \
} while (0)

#define free_3array(a) do {                                    \
    if (a != NULL) {                                           \
        size_t make_3array_loop_counter;                       \
        for (make_3array_loop_counter = 0;                     \
                (a)[make_3array_loop_counter] != NULL;         \
                make_3array_loop_counter++)                    \
            free_matrix((a)[make_3array_loop_counter]);        \
        free_vector(a);                                        \
        a = NULL;                                              \
    }                                                          \
} while (0)

#define make_4array(a, k, l, m, n) do {                             \
    size_t make_4array_loop_counter;                                \
    make_vector(a, (k) + 1);                                        \
    for (make_4array_loop_counter = 0;                              \
            make_4array_loop_counter < (k);                         \
            make_4array_loop_counter++)                             \
        make_3array((a)[make_4array_loop_counter], (l), (m), (n));  \
    (a)[k] = NULL;                                                  \
} while (0)

#define free_4array(a) do {                                    \
    if (a != NULL) {                                           \
        size_t make_4array_loop_counter;                       \
        for (make_4array_loop_counter = 0;                     \
                (a)[make_4array_loop_counter] != NULL;         \
                make_4array_loop_counter++)                    \
            free_3array((a)[make_4array_loop_counter]);        \
        free_vector(a);                                        \
        a = NULL;                                              \
    }                                                          \
} while (0)

#define print_vector(fmt, v, n) do {                           \
    size_t print_vector_loop_counter;                          \
    for (print_vector_loop_counter = 0;                        \
            print_vector_loop_counter < (n);                   \
            print_vector_loop_counter++)                       \
        printf(fmt, (v)[print_vector_loop_counter]);           \
    putchar('\n');                                             \
} while (0)

#define print_matrix(fmt, a, m, n) do {                        \
    size_t print_matrix_loop_counter;                          \
    for (print_matrix_loop_counter = 0;                        \
            print_matrix_loop_counter < (m);                   \
            print_matrix_loop_counter++)                       \
    print_vector(fmt, (a)[print_matrix_loop_counter], (n));    \
} while (0)

#endif /*H_ARRAY_H*/
