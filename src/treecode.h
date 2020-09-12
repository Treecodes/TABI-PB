#ifndef H_TABIPB_TREECODE_STRUCT_H
#define H_TABIPB_TREECODE_STRUCT_H

#include "particles.h"
#include "clusters.h"
#include "interaction_list.h"

class Treecode
{
private:
    class Particles& particles_;
    class Clusters& clusters_;
    const class Tree& tree_;
    const class InteractionList& interaction_list_;
    const class Molecule& molecule_;
    const struct Params& params_;
    
    std::vector<double> potential_;
    
    long int num_iter_;
    
    int gmres_(long int n, const double* b, double* x, long int restrt,
               double* work, long int ldw, double *h, long int ldh,
               long int& iter, double& residual);
    
    void matrix_vector(double alpha, const double* __restrict__ potential_old,
                       double beta,        double* __restrict__ potential_new);
                       
    void precondition(double* z, double* r);
    
    void particle_particle_interact(double* __restrict__ potential,
                              const double* __restrict__ potential_old,
            std::array<std::size_t, 2> target_node_particle_idxs,
            std::array<std::size_t, 2> source_node_particle_idxs);
    
    void particle_cluster_interact(double* __restrict__ potential,
            std::array<std::size_t, 2> target_node_particle_idxs, std::size_t source_node_idx);
                                   
    void cluster_particle_interact(double* __restrict__ potential,
            std::size_t target_node_idx, std::array<std::size_t, 2> source_node_particle_idxs);
            
    void cluster_cluster_interact(double* __restrict__ potential,
            std::size_t target_node_idx, std::size_t source_node_idx);
    
public:
    Treecode(class Particles& particles, class Clusters& clusters,
             const class Tree& tree, const class InteractionList& interaction_list,
             const class Molecule& molecule, const struct Params& params);
    ~Treecode() = default;
    
    void run_GMRES();
    void output();
    
    void copyin_potential_to_device() const;
    void update_potential_on_device() const;
    void update_potential_on_host() const;
    void copyout_potential_to_host() const;
};

#endif /* H_TABIPB_TREECODE_STRUCT_H */
