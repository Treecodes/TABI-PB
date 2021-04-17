#ifndef H_TABIPB_TREE_COMPUTE_STRUCT_H
#define H_TABIPB_TREE_COMPUTE_STRUCT_H

//#include "particles.h"
//#include "interp_pts.h"
#include "tree.h"
#include "interaction_list.h"


class TreeCompute
{
protected:
    const class Tree& source_tree_;
    const class Tree& target_tree_;
    const class InteractionList& interaction_list_;
    
    virtual void particle_particle_interact(std::array<std::size_t, 2> target_node_particle_idxs,
                                            std::array<std::size_t, 2> source_node_particle_idxs) = 0;
    
    virtual void particle_cluster_interact(std::array<std::size_t, 2> target_node_particle_idxs,
                                           std::size_t source_node_idx) = 0;
                                   
    virtual void cluster_particle_interact(std::size_t target_node_idx,
                                           std::array<std::size_t, 2> source_node_particle_idxs) = 0;
            
    virtual void cluster_cluster_interact(std::size_t target_node_idx, std::size_t source_node_idx) = 0;
            
    virtual void upward_pass() = 0;
    virtual void downward_pass() = 0;
    
    virtual void copyin_clusters_to_device() const = 0;
    virtual void delete_clusters_from_device() const = 0;

    
public:
    TreeCompute(const class Tree& source_tree, const class Tree& target_tree,
                const class InteractionList& interaction_list)
        : source_tree_(source_tree), target_tree_(target_tree), interaction_list_(interaction_list) {};
        
    TreeCompute(const class Tree& tree,
                const class InteractionList& interaction_list)
        : source_tree_(tree), target_tree_(tree), interaction_list_(interaction_list) {};
        
    ~TreeCompute() = default;
    
    void run() {
        upward_pass();
#ifdef OPENMP_ENABLED
        #pragma omp parallel for
#endif
        for (std::size_t target_node_idx = 0; target_node_idx < target_tree_.num_nodes(); ++target_node_idx) {
        
            for (auto source_node_idx : interaction_list_.particle_particle(target_node_idx))
                particle_particle_interact(target_tree_.node_particle_idxs(target_node_idx),
                                           source_tree_.node_particle_idxs(source_node_idx));

            for (auto source_node_idx : interaction_list_.particle_cluster(target_node_idx))
                particle_cluster_interact(target_tree_.node_particle_idxs(target_node_idx), source_node_idx);
            
            for (auto source_node_idx : interaction_list_.cluster_particle(target_node_idx))
                cluster_particle_interact(target_node_idx, source_tree_.node_particle_idxs(source_node_idx));
            
            for (auto source_node_idx : interaction_list_.cluster_cluster(target_node_idx))
                cluster_cluster_interact(target_node_idx, source_node_idx);
        }
        downward_pass();
    }
};


#endif /* H_TABIPB_TREE_COMPUTE_STRUCT_H */
