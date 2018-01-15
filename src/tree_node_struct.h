/*
 * C header file for tree node struct
 *
 * C version authored by:
 * Jiahui Chen, Southern Methodist University, Dallas, TX
 * Leighton Wilson, University of Michigan, Ann Arbor, MI
 *
 * Based on package originally written in FORTRAN by:
 * Weihua Geng, Southern Methodist University, Dallas, TX
 * Robery Krasny, University of Michigan, Ann Arbor, MI
 *
 * Last modified by Leighton Wilson, 01/14/2018
 * Created by Leighton Wilson, 01/11/2018
 */

#ifndef H_TREE_NODE_STRUCT_H
#define H_TREE_NODE_STRUCT_H


typedef struct sTreeNode {

    int node_idx;
    int numpar, ibeg, iend;
    double x_min, y_min, z_min;
    double x_max, y_max, z_max;
    double x_mid, y_mid, z_mid;
    double radius, aspect;
    int level, num_children, exist_ms;
    double **ms;
    struct sTreeNode **child;

} TreeNode;

#endif /* H_TREE_NODE_STRUCT */
