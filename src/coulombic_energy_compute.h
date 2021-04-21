#ifndef H_TABIPB_COULOMBIC_ENERGY_COMPUTE_STRUCT_H
#define H_TABIPB_COULOMBIC_ENERGY_COMPUTE_STRUCT_H

#include "timer.h"
#include "molecule.h"
#include "interp_pts.h"
#include "tree_compute.h"

//struct Timers;

class CoulombicEnergyCompute : public TreeCompute
{
private:
    
    const class Molecule& molecule_;
    const class InterpolationPoints& mol_interp_pts_;
    
    const double eps_solute_;
    
    /* Target and Source clusters */
    
    int num_mol_interp_pts_per_node_;
    int num_mol_interp_charges_per_node_;
    std::size_t num_mol_charges_;
    
    int num_mol_interp_potentials_per_node_;
    std::size_t num_mol_potentials_;
    
    std::vector<double> mol_interp_charge_;
    std::vector<double> mol_interp_potential_;
    
    
    /* Coulombic energy */
   
    std::vector<double> coul_eng_vec_; 
    double coulombic_energy_;
    
    
    void particle_particle_interact(std::array<std::size_t, 2> target_node_particle_idxs,
                                    std::array<std::size_t, 2> source_node_particle_idxs) override;
    
    void particle_cluster_interact(std::array<std::size_t, 2> target_node_particle_idxs,
                                   std::size_t source_node_idx) override;
                                   
    void cluster_particle_interact(std::size_t target_node_idx,
                                   std::array<std::size_t, 2> source_node_particle_idxs) override;
            
    void cluster_cluster_interact(std::size_t target_node_idx, std::size_t source_node_idx) override;
            
    void upward_pass() override;
    void downward_pass() override;

    void copyin_clusters_to_device() const override;
    void delete_clusters_from_device() const override;

    
public:

    CoulombicEnergyCompute(const class Molecule& molecule, const class InterpolationPoints& mol_interp_pts,
                      const class Tree& mol_tree, const class InteractionList& interaction_list,
                      double phys_eps_solute);

    ~CoulombicEnergyCompute() = default;
    
    double compute();
    
};

#endif /* H_TABIPB_COULOMBIC_ENERGY_COMPUTE_STRUCT_H */
