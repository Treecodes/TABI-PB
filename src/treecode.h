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
    
    int gmres_(long int n, const double* b, double* x, long int* restrt,
            double* work, long int ldw, double* h, long int ldh,
            long int* iter, double* residual, long int* info);
    
    void matrix_vector(double alpha, double* potential_old, double beta,  double* potential_new);
                       
    void precondition(double* z, double* r);
    
    void particle_particle_interact(double* potential, double* potential_old,
            std::array<std::size_t, 2> target_node_particle_idxs,
            std::array<std::size_t, 2> source_node_particle_idxs);
    
    void particle_cluster_interact(double* potential,
            std::array<std::size_t, 2> target_node_particle_idxs, std::size_t source_node_idx);
                                   
    void cluster_particle_interact(double* potential,
            std::size_t target_node_idx, std::array<std::size_t, 2> source_node_particle_idxs);
            
    void cluster_cluster_interact(double* potential,
            std::size_t target_node_idx, std::size_t source_node_idx);
    
public:
    Treecode(class Particles& particles, class Clusters& clusters,
             const class Tree& tree, const class InteractionList& interaction_list,
             const class Molecule& molecule, const struct Params& params)
        : particles_(particles), clusters_(clusters), tree_(tree),
          interaction_list_(interaction_list), molecule_(molecule), params_(params)
          { potential_.resize(2 * particles_.num()); };
    ~Treecode() = default;
    
    long int run_GMRES();
    void output();
};

#endif /* H_TABIPB_TREECODE_STRUCT_H */
