#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <mpi.h>

#include "../utilities/array.h"
#include "../utilities/tools.h"
#include "../run_params/struct_run_params.h"
#include "../particles/struct_particles.h"

#include "struct_tree_linked_list_node.h"
#include "tree_linked_list.h"

#include "struct_tree.h"
#include "tree.h"


void Tree_Construct(struct Tree **tree_addr, struct Particles *sources, struct RunParams *run_params)
{
    struct TreeLinkedListNode *tree_linked_list = NULL;
    
    double xyzminmax[6];
    int numnodes = 0;
    int numleaves = 0;
    int max_depth = 1;
    
    int min_leaf_size = INT_MAX;
    int max_leaf_size = 0;
    
    xyzminmax[0] = minval(sources->x, sources->num);
    xyzminmax[1] = maxval(sources->x, sources->num);
    xyzminmax[2] = minval(sources->y, sources->num);
    xyzminmax[3] = maxval(sources->y, sources->num);
    xyzminmax[4] = minval(sources->z, sources->num);
    xyzminmax[5] = maxval(sources->z, sources->num);
    
    TreeLinkedList_Construct(&tree_linked_list, NULL, sources, 1, sources->num,
                    run_params->max_per_source_leaf, xyzminmax, &numnodes, &numleaves,
                    &min_leaf_size, &max_leaf_size, &max_depth, 0);
    
    TreeLinkedList_SetIndex(tree_linked_list, 0);
    


    Tree_Alloc(tree_addr, numnodes);
    Tree_Fill(*tree_addr, tree_linked_list);
    (*tree_addr)->numleaves = numleaves;
    
    (*tree_addr)->min_leaf_size = min_leaf_size;
    (*tree_addr)->max_leaf_size = max_leaf_size;
    (*tree_addr)->max_depth = max_depth;
    
    Tree_Set_Leaves_and_Levels(*tree_addr);

    TreeLinkedList_Free(&tree_linked_list);

    return;
}


void Tree_Set_Leaves_and_Levels(struct Tree *tree)
{

    /* Creating levels list for the downpass */
    make_matrix(tree->levels_list, tree->max_depth, 20);
    make_vector(tree->levels_list_num, tree->max_depth);
    for (int i = 0; i < tree->max_depth; ++i) tree->levels_list_num[i]=0;

    make_vector(tree->leaves_list, 50);
    tree->leaves_list_num = 0;

    int *sizeof_levels_list = NULL;
    make_vector(sizeof_levels_list, tree->max_depth);
    for (int i = 0; i < tree->max_depth; ++i) sizeof_levels_list[i]=20;

    int sizeof_leaves_list = 50;

    Tree_Fill_Levels(tree, 0, 0, sizeof_levels_list, &sizeof_leaves_list);
    free_vector(sizeof_levels_list);

    return;
}


void Tree_Fill_Levels(struct Tree *tree, int idx, int level, int *sizeof_levels_list, int *sizeof_leaves_list)
{

    if (tree->num_children[idx] == 0) {
        if (tree->leaves_list_num >= *sizeof_leaves_list) {
            *sizeof_leaves_list *= 1.5;
            tree->leaves_list = realloc_vector(tree->leaves_list, *sizeof_leaves_list);
        }

        tree->leaves_list[tree->leaves_list_num] = idx;
        tree->leaves_list_num++;

    } else {

        if (tree->levels_list_num[level] >= sizeof_levels_list[level]) {
            sizeof_levels_list[level] *= 1.5;
            tree->levels_list[level] = realloc_vector(tree->levels_list[level], sizeof_levels_list[level]);
        }

        tree->levels_list[level][tree->levels_list_num[level]] = idx;
        tree->levels_list_num[level]++;

        for (int i = 0; i < tree->num_children[idx]; i++)
            Tree_Fill_Levels(tree, tree->children[8*idx + i], level+1, sizeof_levels_list, sizeof_leaves_list);
    }


    return;
}


void Tree_Alloc(struct Tree **tree_addr, int length)
{
    *tree_addr = malloc(sizeof(struct Tree));
    struct Tree *tree = *tree_addr;

    tree->numnodes = length;
    make_vector(tree->ibeg, length);
    make_vector(tree->iend, length);
    make_vector(tree->numpar, length);
    make_vector(tree->x_mid, length);
    make_vector(tree->y_mid, length);
    make_vector(tree->z_mid, length);
    make_vector(tree->x_min, length);
    make_vector(tree->y_min, length);
    make_vector(tree->z_min, length);
    make_vector(tree->x_max, length);
    make_vector(tree->y_max, length);
    make_vector(tree->z_max, length);
    make_vector(tree->cluster_ind, length);
    make_vector(tree->radius, length);

    make_vector(tree->num_children, length);
    make_vector(tree->children, 8*length);
    make_vector(tree->parent, length);

    tree->levels_list = NULL;
    tree->levels_list_num = NULL;
    tree->leaves_list = NULL;

    
    return;
}   /* END of function allocate_tree */



void Tree_Free(struct Tree **tree_addr)
{
    struct Tree *tree = *tree_addr;
    
    if (tree != NULL) {
        free_vector(tree->ibeg);
        free_vector(tree->iend);
        free_vector(tree->numpar);
        free_vector(tree->x_mid);
        free_vector(tree->y_mid);
        free_vector(tree->z_mid);
        free_vector(tree->x_min);
        free_vector(tree->y_min);
        free_vector(tree->z_min);
        free_vector(tree->x_max);
        free_vector(tree->y_max);
        free_vector(tree->z_max);
        free_vector(tree->cluster_ind);
        free_vector(tree->radius);

        free_vector(tree->num_children);
        free_vector(tree->children);
        free_vector(tree->parent);

        if (tree->levels_list != NULL) free_matrix(tree->levels_list);
        if (tree->levels_list_num != NULL) free_vector(tree->levels_list_num);
        if (tree->leaves_list != NULL) free_vector(tree->leaves_list);

        free(tree);
    }

    tree = NULL;

    return;
}   /* END of function allocate_tree */



void Tree_Fill(struct Tree *tree, struct TreeLinkedListNode *p)
{
    tree->x_mid[p->node_index] = p->x_mid;
    tree->y_mid[p->node_index] = p->y_mid;
    tree->z_mid[p->node_index] = p->z_mid;
    
    tree->x_min[p->node_index] = p->x_min;
    tree->y_min[p->node_index] = p->y_min;
    tree->z_min[p->node_index] = p->z_min;
    
    tree->x_max[p->node_index] = p->x_max;
    tree->y_max[p->node_index] = p->y_max;
    tree->z_max[p->node_index] = p->z_max;
    
    tree->ibeg[p->node_index] = p->ibeg;
    tree->iend[p->node_index] = p->iend;
    tree->numpar[p->node_index] = p->numpar;
    tree->radius[p->node_index] = p->radius;
    tree->cluster_ind[p->node_index] = p->node_index;

    tree->num_children[p->node_index] = p->num_children;

    if (p->parent != NULL)
        tree->parent[p->node_index] = (p->parent)->node_index;
    else
        tree->parent[p->node_index] = -1;

    for (int i = 0; i < p->num_children; i++) {
        tree->children[8*p->node_index+i] = (p->child[i])->node_index;
        Tree_Fill(tree, p->child[i]);
    }
    
    return;
} /* END of function Tree_CreateArray */
