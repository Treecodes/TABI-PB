#ifndef H_TABIPB_TREECODE_STRUCT_H
#define H_TABIPB_TREECODE_STRUCT_H

#include "timer.h"
#include "output.h"
#include "elements.h"
#include "interp_pts.h"
#include "interaction_list.h"

struct Timers_BoundaryElement;

class BoundaryElement
{
private:
    class Elements& elements_;
    const class InterpolationPoints& interp_pts_;
    const class Tree& tree_;
    const class InteractionList& interaction_list_;
    const class Molecule& molecule_;
    const struct Params& params_;
    class Output& output_;
    struct Timers_BoundaryElement& timers_;
    
    std::vector<double> potential_;
    
    /* cluster specific data */
    int num_charges_per_node_;
    std::size_t num_charges_;
    
    std::vector<double> interp_charge_;
    std::vector<double> interp_charge_dx_;
    std::vector<double> interp_charge_dy_;
    std::vector<double> interp_charge_dz_;
    
    std::vector<double> interp_potential_;
    std::vector<double> interp_potential_dx_;
    std::vector<double> interp_potential_dy_;
    std::vector<double> interp_potential_dz_;
    
    /* output */
    double solvation_energy_;
    double free_energy_;
    double coulombic_energy_;
    
    double pot_min_;
    double pot_max_;
    double pot_normal_min_;
    double pot_normal_max_;
    
    int gmres_(long int n, const double* b, double* x, long int restrt,
               double* work, long int ldw, double *h, long int ldh,
               long int& iter, double& residual);
    
    void matrix_vector(double alpha, const double* __restrict potential_old,
                       double beta,        double* __restrict potential_new);
                       
    void precondition_diagonal(double* z, double* r);
    void precondition_block(double* z, double* r);
    
    void particle_particle_interact(double* __restrict potential,
                              const double* __restrict potential_old,
            std::array<std::size_t, 2> target_node_particle_idxs,
            std::array<std::size_t, 2> source_node_particle_idxs);
    
    void particle_cluster_interact(double* __restrict potential,
            std::array<std::size_t, 2> target_node_particle_idxs, std::size_t source_node_idx);
                                   
    void cluster_particle_interact(double* __restrict potential,
            std::size_t target_node_idx, std::array<std::size_t, 2> source_node_particle_idxs);
            
    void cluster_cluster_interact(double* __restrict potential,
            std::size_t target_node_idx, std::size_t source_node_idx);
            
    void upward_pass();
    void downward_pass(double* __restrict potential);
    
    void clear_cluster_charges();
    void clear_cluster_potentials();
    void copyin_clusters_to_device() const;
    void delete_clusters_from_device() const;
    
    const std::array<std::size_t, 2> cluster_charges_idxs(std::size_t node_idx) const {
        return std::array<std::size_t, 2> {num_charges_per_node_ *  node_idx,
                                           num_charges_per_node_ * (node_idx + 1)};
    };

    
public:
    BoundaryElement(class Elements& elements, const class InterpolationPoints& interp_pts,
             const class Tree& tree, const class InteractionList& interaction_list,
             const class Molecule& molecule, const struct Params& params, class Output& output,
             struct Timers_BoundaryElement& timers);
    ~BoundaryElement() = default;
    
    void run_GMRES();
    //void finalize();

};


struct Timers_BoundaryElement
{
    Timer ctor;
    Timer run_GMRES;
    
    Timer clear_charges;
    Timer clear_potentials;

    Timer matrix_vector;
    Timer precondition;
    
    Timer particle_particle_interact;
    Timer particle_cluster_interact;
    Timer cluster_particle_interact;
    Timer cluster_cluster_interact;
    Timer upward_pass;
    Timer downward_pass;
    
    Timer clear_cluster_charges;
    Timer clear_cluster_potentials;
    Timer copyin_clusters_to_device;
    Timer delete_clusters_from_device;

    void print() const;
    std::string get_durations() const;
    std::string get_headers() const;

    Timers_BoundaryElement() = default;
    ~Timers_BoundaryElement() = default;
};

#endif /* H_TABIPB_TREECODE_STRUCT_H */
