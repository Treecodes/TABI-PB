#include <cmath>
#include <algorithm>
#include <vector>

#include "constants.h"
#include "coulombic_energy_compute.h"


CoulombicEnergyCompute::CoulombicEnergyCompute(const class Molecule& molecule,
                      const class InterpolationPoints& mol_interp_pts, const class Tree& mol_tree,
                      const class InteractionList& interaction_list, double phys_eps_solute)
    : TreeCompute(mol_tree, interaction_list),
      molecule_(molecule), mol_interp_pts_(mol_interp_pts), eps_solute_(phys_eps_solute)
{
//    timers_.ctor.start();
    
    /* Target and Source clusters */
    
    num_mol_interp_pts_per_node_        = mol_interp_pts_.num_interp_pts_per_node();
    num_mol_interp_charges_per_node_    = std::pow(num_mol_interp_pts_per_node_, 3);
    num_mol_charges_                    = source_tree_.num_nodes() * num_mol_interp_charges_per_node_;
    
    num_mol_interp_potentials_per_node_ = num_mol_interp_charges_per_node_;
    num_mol_potentials_                 = num_mol_charges_;
    
    mol_interp_charge_.assign(num_mol_charges_, 0.);
    mol_interp_potential_.assign(num_mol_potentials_, 0.);

    /* Coulombic energy */

    coul_eng_vec_.resize(1);
    coul_eng_vec_[0] = 0.;

    coulombic_energy_ = 0.;

//    timers_.ctor.stop();
}


double CoulombicEnergyCompute::compute()
{
    CoulombicEnergyCompute::copyin_clusters_to_device();
    CoulombicEnergyCompute::run();
    CoulombicEnergyCompute::delete_clusters_from_device();
    
    coulombic_energy_ = coul_eng_vec_[0] / 2;
    
    return coulombic_energy_;
}


void CoulombicEnergyCompute::particle_particle_interact(std::array<std::size_t, 2> target_node_idxs,
                                                        std::array<std::size_t, 2> source_node_idxs)
{
//    timers_.particle_particle_interact.start();

    /* Targets */
    
    std::size_t target_node_begin      = target_node_idxs[0];
    std::size_t target_node_end        = target_node_idxs[1];

    /* Sources */
    
    std::size_t source_node_begin      = source_node_idxs[0];
    std::size_t source_node_end        = source_node_idxs[1];
    
    const double* __restrict mol_x_ptr = molecule_.x_ptr();
    const double* __restrict mol_y_ptr = molecule_.y_ptr();
    const double* __restrict mol_z_ptr = molecule_.z_ptr();

    const double* __restrict mol_q_ptr = molecule_.charge_ptr();

    /* Potential */

    double* __restrict coul_eng_ptr = coul_eng_vec_.data();

#ifdef OPENACC_ENABLED
    int stream_id = std::rand() % 3;
    #pragma acc parallel loop async(stream_id) present(mol_x_ptr, mol_y_ptr, mol_z_ptr, mol_q_ptr, \
                                                       coul_eng_ptr)
#endif
    for (std::size_t j = target_node_begin; j < target_node_end; ++j) {
        
        double target_x = mol_x_ptr[j];
        double target_y = mol_y_ptr[j];
        double target_z = mol_z_ptr[j];
        double target_q = mol_q_ptr[j];
        
        double pot_temp = 0.;
        
#ifdef OPENACC_ENABLED
        #pragma acc loop reduction(+:pot_temp)
#endif
        for (std::size_t k = source_node_begin; k < source_node_end; ++k) {

            double dx = target_x - mol_x_ptr[k];
            double dy = target_y - mol_y_ptr[k];
            double dz = target_z - mol_z_ptr[k];
            double r  = dx*dx + dy*dy + dz*dz;

            if (r > 0) pot_temp += target_q * mol_q_ptr[k] / eps_solute_ / std::sqrt(r);
        }

#ifdef OPENACC_ENABLED
        #pragma acc atomic update
#elif  OPENMP_ENABLED
        #pragma omp atomic update
#endif
        coul_eng_ptr[0] += pot_temp;
    }

//    timers_.particle_particle_interact.stop();
}


void CoulombicEnergyCompute::particle_cluster_interact(std::array<std::size_t, 2> target_node_idxs,
                                                       std::size_t source_node_idx)
{
//    timers_.particle_cluster_interact.start();
    
    /* Targets */
    
    std::size_t target_node_begin          = target_node_idxs[0];
    std::size_t target_node_end            = target_node_idxs[1];
    
    const double* __restrict mol_x_ptr = molecule_.x_ptr();
    const double* __restrict mol_y_ptr = molecule_.y_ptr();
    const double* __restrict mol_z_ptr = molecule_.z_ptr();

    const double* __restrict mol_q_ptr = molecule_.charge_ptr();

    /* Sources */
    
    int num_mol_interp_pts_per_node                 = num_mol_interp_pts_per_node_;
    int num_mol_interp_charges_per_node             = num_mol_interp_charges_per_node_;

    std::size_t source_cluster_interp_pts_begin     = source_node_idx * num_mol_interp_pts_per_node_;
    std::size_t source_cluster_interp_charges_begin = source_node_idx * num_mol_interp_charges_per_node_;
    
    const double* __restrict mol_clusters_x_ptr     = mol_interp_pts_.interp_x_ptr();
    const double* __restrict mol_clusters_y_ptr     = mol_interp_pts_.interp_y_ptr();
    const double* __restrict mol_clusters_z_ptr     = mol_interp_pts_.interp_z_ptr();
    
    const double* __restrict mol_clusters_q_ptr     = mol_interp_charge_.data();

    /* Potential */

    double* __restrict coul_eng_ptr = coul_eng_vec_.data();
    
    
#ifdef OPENACC_ENABLED
    int stream_id = std::rand() % 3;
    #pragma acc parallel loop async(stream_id) present(mol_x_ptr, mol_y_ptr, mol_z_ptr, mol_q_ptr, \
                                      mol_clusters_x_ptr, mol_clusters_y_ptr, mol_clusters_z_ptr, \
                                      mol_clusters_q_ptr, coul_eng_ptr)
#endif
    for (std::size_t j = target_node_begin; j < target_node_end; ++j) {

        double target_x = mol_x_ptr[j];
        double target_y = mol_y_ptr[j];
        double target_z = mol_z_ptr[j];
        double target_q = mol_q_ptr[j];
        
        double pot_temp = 0.;
        
#ifdef OPENACC_ENABLED
        #pragma acc loop collapse(3) reduction(+:pot_temp)
#endif
        for (int k1 = 0; k1 < num_mol_interp_pts_per_node; ++k1) {
        for (int k2 = 0; k2 < num_mol_interp_pts_per_node; ++k2) {
        for (int k3 = 0; k3 < num_mol_interp_pts_per_node; ++k3) {
                
            std::size_t kk = source_cluster_interp_charges_begin
                           + k1 * num_mol_interp_pts_per_node * num_mol_interp_pts_per_node
                           + k2 * num_mol_interp_pts_per_node + k3;

            double dx = target_x - mol_clusters_x_ptr[source_cluster_interp_pts_begin + k1];
            double dy = target_y - mol_clusters_y_ptr[source_cluster_interp_pts_begin + k2];
            double dz = target_z - mol_clusters_z_ptr[source_cluster_interp_pts_begin + k3];

            pot_temp += mol_clusters_q_ptr[kk] / eps_solute_ / std::sqrt(dx*dx + dy*dy + dz*dz);
        }
        }
        }
        
        pot_temp *= target_q;

#ifdef OPENACC_ENABLED
        #pragma acc atomic update
#elif  OPENMP_ENABLED
        #pragma omp atomic update
#endif
        coul_eng_ptr[0] += pot_temp;
    }

//    timers_.particle_cluster_interact.stop();
}


void CoulombicEnergyCompute::cluster_particle_interact(std::size_t target_node_idx,
                                                       std::array<std::size_t, 2> source_node_idxs)
{
//    timers_.cluster_particle_interact.start();

    /* Targets */
    
    int num_mol_interp_pts_per_node                    = num_mol_interp_pts_per_node_;
    int num_mol_interp_potentials_per_node             = num_mol_interp_potentials_per_node_;

    std::size_t target_cluster_interp_pts_begin        = target_node_idx * num_mol_interp_pts_per_node_;
    std::size_t target_cluster_interp_potentials_begin = target_node_idx * num_mol_interp_potentials_per_node_;

    const double* __restrict mol_clusters_x_ptr = mol_interp_pts_.interp_x_ptr();
    const double* __restrict mol_clusters_y_ptr = mol_interp_pts_.interp_y_ptr();
    const double* __restrict mol_clusters_z_ptr = mol_interp_pts_.interp_z_ptr();

    double*       __restrict mol_clusters_p_ptr = mol_interp_potential_.data();


    /* Sources */
    
    std::size_t source_node_begin      = source_node_idxs[0];
    std::size_t source_node_end        = source_node_idxs[1];

    const double* __restrict mol_x_ptr = molecule_.x_ptr();
    const double* __restrict mol_y_ptr = molecule_.y_ptr();
    const double* __restrict mol_z_ptr = molecule_.z_ptr();

    const double* __restrict mol_q_ptr = molecule_.charge_ptr();


#ifdef OPENACC_ENABLED
    int stream_id = std::rand() % 3;
    #pragma acc parallel loop collapse(3) async(stream_id) present(mol_x_ptr, mol_y_ptr, mol_z_ptr, mol_q_ptr, \
                    mol_clusters_x_ptr, mol_clusters_y_ptr, mol_clusters_z_ptr, mol_clusters_p_ptr)
#endif
    for (int j1 = 0; j1 < num_mol_interp_pts_per_node; ++j1) {
    for (int j2 = 0; j2 < num_mol_interp_pts_per_node; ++j2) {
    for (int j3 = 0; j3 < num_mol_interp_pts_per_node; ++j3) {
    
        std::size_t jj = target_cluster_interp_potentials_begin
                       + j1 * num_mol_interp_pts_per_node * num_mol_interp_pts_per_node
                       + j2 * num_mol_interp_pts_per_node + j3;

        double target_x = mol_clusters_x_ptr[target_cluster_interp_pts_begin + j1];
        double target_y = mol_clusters_y_ptr[target_cluster_interp_pts_begin + j2];
        double target_z = mol_clusters_z_ptr[target_cluster_interp_pts_begin + j3];
        
        double pot_temp = 0.;
    
#ifdef OPENACC_ENABLED
        #pragma acc loop reduction(+:pot_temp)
#endif
        for (std::size_t k = source_node_begin; k < source_node_end; ++k) {

            double dx = target_x - mol_x_ptr[k];
            double dy = target_y - mol_y_ptr[k];
            double dz = target_z - mol_z_ptr[k];

            pot_temp += mol_q_ptr[k] / eps_solute_ / std::sqrt(dx*dx + dy*dy + dz*dz);;
        }
    
#ifdef OPENACC_ENABLED
        #pragma acc atomic update
#elif  OPENMP_ENABLED
        #pragma omp atomic update
#endif
        mol_clusters_p_ptr[jj] += pot_temp;

    }
    }
    }

//    timers_.cluster_particle_interact.stop();
}


void CoulombicEnergyCompute::cluster_cluster_interact(std::size_t target_node_idx,
                                                      std::size_t source_node_idx)
{
//    timers_.cluster_cluster_interact.start();
        
    /* Targets */
    
    int num_mol_interp_potentials_per_node            = num_mol_interp_potentials_per_node_;
    
    std::size_t target_cluster_interp_pts_begin        = target_node_idx * num_mol_interp_pts_per_node_;
    std::size_t target_cluster_interp_potentials_begin = target_node_idx * num_mol_interp_potentials_per_node_;

    double*       __restrict mol_clusters_p_ptr       = mol_interp_potential_.data();
    
    
    /* Sources */
    
    int num_mol_interp_pts_per_node                 = num_mol_interp_pts_per_node_;
    int num_mol_interp_charges_per_node             = num_mol_interp_charges_per_node_;

    std::size_t source_cluster_interp_pts_begin     = source_node_idx * num_mol_interp_pts_per_node;
    std::size_t source_cluster_interp_charges_begin = source_node_idx * num_mol_interp_charges_per_node;
    
    const double* __restrict mol_clusters_x_ptr     = mol_interp_pts_.interp_x_ptr();
    const double* __restrict mol_clusters_y_ptr     = mol_interp_pts_.interp_y_ptr();
    const double* __restrict mol_clusters_z_ptr     = mol_interp_pts_.interp_z_ptr();
    
    const double* __restrict mol_clusters_q_ptr     = mol_interp_charge_.data();


#ifdef OPENACC_ENABLED
    int stream_id = std::rand() % 3;
    #pragma acc parallel loop collapse(3) async(stream_id) present(mol_clusters_x_ptr, mol_clusters_y_ptr, \
                    mol_clusters_z_ptr, mol_clusters_q_ptr,  mol_clusters_p_ptr)
#endif
    for (int j1 = 0; j1 < num_mol_interp_pts_per_node; j1++) {
    for (int j2 = 0; j2 < num_mol_interp_pts_per_node; j2++) {
    for (int j3 = 0; j3 < num_mol_interp_pts_per_node; j3++) {
    
        std::size_t jj = target_cluster_interp_potentials_begin
                       + j1 * num_mol_interp_pts_per_node * num_mol_interp_pts_per_node
                       + j2 * num_mol_interp_pts_per_node + j3;

        double target_x = mol_clusters_x_ptr[target_cluster_interp_pts_begin + j1];
        double target_y = mol_clusters_y_ptr[target_cluster_interp_pts_begin + j2];
        double target_z = mol_clusters_z_ptr[target_cluster_interp_pts_begin + j3];
        
        double pot_temp = 0.;
    
#ifdef OPENACC_ENABLED
        #pragma acc loop collapse(3) reduction(+:pot_temp)
#endif
        for (int k1 = 0; k1 < num_mol_interp_pts_per_node; k1++) {
        for (int k2 = 0; k2 < num_mol_interp_pts_per_node; k2++) {
        for (int k3 = 0; k3 < num_mol_interp_pts_per_node; k3++) {
            
            std::size_t kk = source_cluster_interp_charges_begin
                           + k1 * num_mol_interp_pts_per_node * num_mol_interp_pts_per_node
                           + k2 * num_mol_interp_pts_per_node + k3;

            double dx = target_x - mol_clusters_x_ptr[source_cluster_interp_pts_begin + k1];
            double dy = target_y - mol_clusters_y_ptr[source_cluster_interp_pts_begin + k2];
            double dz = target_z - mol_clusters_z_ptr[source_cluster_interp_pts_begin + k3];

            pot_temp += mol_clusters_q_ptr[kk] / eps_solute_ / std::sqrt(dx*dx + dy*dy + dz*dz);;
        }
        }
        }
    
#ifdef OPENACC_ENABLED
        #pragma acc atomic update
#elif  OPENMP_ENABLED
        #pragma omp atomic update
#endif
        mol_clusters_p_ptr[jj] += pot_temp;
    }
    }
    }

//    timers_.cluster_cluster_interact.stop();
}


void CoulombicEnergyCompute::upward_pass()
{
//    timers_.upward_pass.start();

    int num_mol_interp_pts_per_node     = num_mol_interp_pts_per_node_;
    int num_mol_interp_charges_per_node = num_mol_interp_charges_per_node_;
    
    const double* __restrict mol_x_ptr = molecule_.x_ptr();
    const double* __restrict mol_y_ptr = molecule_.y_ptr();
    const double* __restrict mol_z_ptr = molecule_.z_ptr();
    
    const double* __restrict mol_q_ptr = molecule_.charge_ptr();
    
    const double* __restrict mol_clusters_x_ptr = mol_interp_pts_.interp_x_ptr();
    const double* __restrict mol_clusters_y_ptr = mol_interp_pts_.interp_y_ptr();
    const double* __restrict mol_clusters_z_ptr = mol_interp_pts_.interp_z_ptr();
    
    double*       __restrict mol_clusters_q_ptr = mol_interp_charge_.data();
    
        
    std::vector<double> weights (num_mol_interp_pts_per_node);
    double* weights_ptr = weights.data();
    int weights_num = weights.size();
    
    for (int i = 0; i < weights_num; ++i) {
        weights[i] = ((i % 2 == 0)? 1 : -1);
        if (i == 0 || i == weights_num-1) weights[i] = ((i % 2 == 0)? 1 : -1) * 0.5;
    }

    
#ifdef OPENACC_ENABLED
    #pragma acc enter data copyin(weights_ptr[0:weights_num])
#endif
    
    for (std::size_t node_idx = 0; node_idx < source_tree_.num_nodes(); ++node_idx) {
        
        auto particle_idxs = source_tree_.node_particle_idxs(node_idx);
        
        std::size_t node_interp_pts_start = node_idx * num_mol_interp_pts_per_node;
        std::size_t node_charges_start    = node_idx * num_mol_interp_charges_per_node;
        
        std::size_t particle_start = particle_idxs[0];
        std::size_t num_particles  = particle_idxs[1] - particle_idxs[0];
        
        std::vector<int> exact_idx_x(num_particles);
        std::vector<int> exact_idx_y(num_particles);
        std::vector<int> exact_idx_z(num_particles);
        std::vector<double> denominator(num_particles);
        
        int* exact_idx_x_ptr = exact_idx_x.data();
        int* exact_idx_y_ptr = exact_idx_y.data();
        int* exact_idx_z_ptr = exact_idx_z.data();
        double* denominator_ptr = denominator.data();
        
#ifdef OPENACC_ENABLED
#pragma acc kernels present(mol_x_ptr, mol_y_ptr, mol_z_ptr, mol_q_ptr, \
                            mol_clusters_x_ptr, mol_clusters_y_ptr, mol_clusters_z_ptr, \
                            mol_clusters_q_ptr, weights_ptr) \
                  create(exact_idx_x_ptr[0:num_particles], exact_idx_y_ptr[0:num_particles], \
                         exact_idx_z_ptr[0:num_particles], denominator_ptr[0:num_particles])
#endif
        {

#ifdef OPENACC_ENABLED
        #pragma acc loop vector(32) independent
#endif
        for (std::size_t i = 0; i < num_particles; ++i) {
            exact_idx_x_ptr[i] = -1;
            exact_idx_y_ptr[i] = -1;
            exact_idx_z_ptr[i] = -1;
        }

#ifdef OPENACC_ENABLED
        #pragma acc loop independent
#endif
        for (std::size_t i = 0; i < num_particles; ++i) {
        
            double denominator_x = 0.;
            double denominator_y = 0.;
            double denominator_z = 0.;
            
            double xx = mol_x_ptr[particle_start + i];
            double yy = mol_y_ptr[particle_start + i];
            double zz = mol_z_ptr[particle_start + i];

            // because there's a reduction over exact_idx[i], this loop carries a
            // backward dependence and won't actually parallelize
#ifdef OPENACC_ENABLED
            #pragma acc loop reduction(+:denominator_x,denominator_y,denominator_z)
#endif
            for (int j = 0; j < num_mol_interp_pts_per_node; ++j) {
            
                double dist_x = xx - mol_clusters_x_ptr[node_interp_pts_start + j];
                double dist_y = yy - mol_clusters_y_ptr[node_interp_pts_start + j];
                double dist_z = zz - mol_clusters_z_ptr[node_interp_pts_start + j];
                
                denominator_x += weights_ptr[j] / dist_x;
                denominator_y += weights_ptr[j] / dist_y;
                denominator_z += weights_ptr[j] / dist_z;
                
                if (std::abs(dist_x) < std::numeric_limits<double>::min()) exact_idx_x_ptr[i] = j;
                if (std::abs(dist_y) < std::numeric_limits<double>::min()) exact_idx_y_ptr[i] = j;
                if (std::abs(dist_z) < std::numeric_limits<double>::min()) exact_idx_z_ptr[i] = j;
            }
            
            denominator_ptr[i] = 1.0;
            if (exact_idx_x_ptr[i] == -1) denominator_ptr[i] /= denominator_x;
            if (exact_idx_y_ptr[i] == -1) denominator_ptr[i] /= denominator_y;
            if (exact_idx_z_ptr[i] == -1) denominator_ptr[i] /= denominator_z;
        }

#ifdef OPENACC_ENABLED
        #pragma acc loop collapse(3) independent
#endif
        for (int k1 = 0; k1 < num_mol_interp_pts_per_node; ++k1) {
        for (int k2 = 0; k2 < num_mol_interp_pts_per_node; ++k2) {
        for (int k3 = 0; k3 < num_mol_interp_pts_per_node; ++k3) {
        
            std::size_t kk = node_charges_start
                   + k1 * num_mol_interp_pts_per_node * num_mol_interp_pts_per_node
                   + k2 * num_mol_interp_pts_per_node + k3;
                   
            double cx = mol_clusters_x_ptr[node_interp_pts_start + k1];
            double w1 = weights_ptr[k1];

            double cy = mol_clusters_y_ptr[node_interp_pts_start + k2];
            double w2 = weights_ptr[k2];
            
            double cz = mol_clusters_z_ptr[node_interp_pts_start + k3];
            double w3 = weights_ptr[k3];
            
            double q_temp = 0.;
            
#ifdef OPENACC_ENABLED
            #pragma acc loop reduction(+:q_temp)
#endif
            for (std::size_t i = 0; i < num_particles; i++) {  // loop over source points
            
                double dist_x = mol_x_ptr[particle_start + i] - cx;
                double dist_y = mol_y_ptr[particle_start + i] - cy;
                double dist_z = mol_z_ptr[particle_start + i] - cz;
                
                double numerator = 1.;

                // If exact_idx[i] == -1, then no issues.
                // If exact_idx[i] != -1, then we want to zero out terms EXCEPT when exactInd=k1.
                if (exact_idx_x_ptr[i] == -1) {
                    numerator *= w1 / dist_x;
                } else {
                    if (exact_idx_x_ptr[i] != k1) numerator *= 0.;
                }

                if (exact_idx_y_ptr[i] == -1) {
                    numerator *= w2 / dist_y;
                } else {
                    if (exact_idx_y_ptr[i] != k2) numerator *= 0.;
                }

                if (exact_idx_z_ptr[i] == -1) {
                    numerator *= w3 / dist_z;
                } else {
                    if (exact_idx_z_ptr[i] != k3) numerator *= 0.;
                }

                q_temp += mol_q_ptr[particle_start + i] * numerator * denominator_ptr[i];
            }
            
            mol_clusters_q_ptr[kk] += q_temp;
        }
        }
        }
        
        } // end parallel region
    } // end loop over nodes
#ifdef OPENACC_ENABLED
    #pragma acc exit data delete(weights_ptr[0:weights_num])
#endif

//    timers_.upward_pass.stop();
}


void CoulombicEnergyCompute::downward_pass()
{
/*
//    timers_.downward_pass.start();

    int num_mol_interp_pts_per_node        = num_mol_interp_pts_per_node_;
    int num_mol_interp_potentials_per_node = num_mol_interp_potentials_per_node_;
    
    const double* __restrict mol_x_ptr = molecule_.x_ptr();
    const double* __restrict mol_y_ptr = molecule_.y_ptr();
    const double* __restrict mol_z_ptr = molecule_.z_ptr();
    
    const double* __restrict mol_q_ptr = molecule_.charge_ptr();
        
    const double* __restrict mol_clusters_x_ptr = mol_interp_pts_.interp_x_ptr();
    const double* __restrict mol_clusters_y_ptr = mol_interp_pts_.interp_y_ptr();
    const double* __restrict mol_clusters_z_ptr = mol_interp_pts_.interp_z_ptr();
    
    const double* __restrict mol_clusters_p_ptr = mol_interp_potential_.data();

    double* __restrict coul_eng_ptr = coul_eng_vec_.data();
    
    std::vector<double> weights (num_mol_interp_pts_per_node);
    double* weights_ptr = weights.data();
    int weights_num = weights.size();
    
    for (int i = 0; i < weights_num; ++i) {
        weights[i] = ((i % 2 == 0)? 1 : -1);
        if (i == 0 || i == weights_num-1) weights[i] = ((i % 2 == 0)? 1 : -1) * 0.5;
    }

    
#ifdef OPENACC_ENABLED
#pragma acc enter data copyin(weights_ptr[0:weights_num])
#endif
    
    for (std::size_t node_idx = 0; node_idx < target_tree_.num_nodes(); ++node_idx) {
        
        auto particle_idxs = target_tree_.node_particle_idxs(node_idx);
        std::size_t node_interp_pts_start = node_idx * num_mol_interp_pts_per_node;
        std::size_t node_potentials_start = node_idx * num_mol_interp_potentials_per_node;
        
        std::size_t particle_start = particle_idxs[0];
        std::size_t num_particles  = particle_idxs[1] - particle_idxs[0];

#ifdef OPENACC_ENABLED
#pragma acc parallel loop present(mol_x_ptr, mol_y_ptr, mol_z_ptr, mol_clusters_p_ptr, \
                                  coul_eng_ptr, weights_ptr)
#endif
        for (std::size_t i = 0; i < num_particles; ++i) {
        
            double denominator_x = 0.;
            double denominator_y = 0.;
            double denominator_z = 0.;
            
            int exact_idx_x = -1;
            int exact_idx_y = -1;
            int exact_idx_z = -1;
            
            double xx = mol_x_ptr[particle_start + i];
            double yy = mol_y_ptr[particle_start + i];
            double zz = mol_z_ptr[particle_start + i];
            double qq = mol_q_ptr[particle_start + i];
            
#ifdef OPENACC_ENABLED
            #pragma acc loop reduction(+:denominator_x,denominator_y,denominator_z) \
                             reduction(max:exact_idx_x,exact_idx_y,exact_idx_z)
#endif
            for (int j = 0; j < num_mol_interp_pts_per_node; ++j) {
            
                double dist_x = xx - mol_clusters_x_ptr[node_interp_pts_start + j];
                double dist_y = yy - mol_clusters_y_ptr[node_interp_pts_start + j];
                double dist_z = zz - mol_clusters_z_ptr[node_interp_pts_start + j];
                
                denominator_x += weights_ptr[j] / dist_x;
                denominator_y += weights_ptr[j] / dist_y;
                denominator_z += weights_ptr[j] / dist_z;
                
                if (std::abs(dist_x) < std::numeric_limits<double>::min()) exact_idx_x = j;
                if (std::abs(dist_y) < std::numeric_limits<double>::min()) exact_idx_y = j;
                if (std::abs(dist_z) < std::numeric_limits<double>::min()) exact_idx_z = j;
            }
            
            double denominator = 1.;
            if (exact_idx_x == -1) denominator /= denominator_x;
            if (exact_idx_y == -1) denominator /= denominator_y;
            if (exact_idx_z == -1) denominator /= denominator_z;

            double pot_temp = 0.;
            
#ifdef OPENACC_ENABLED
            #pragma acc loop collapse(3) reduction(+:pot_temp)
#endif
            for (int k1 = 0; k1 < num_mol_interp_pts_per_node; ++k1) {
            for (int k2 = 0; k2 < num_mol_interp_pts_per_node; ++k2) {
            for (int k3 = 0; k3 < num_mol_interp_pts_per_node; ++k3) {
                    
                std::size_t kk = node_potentials_start
                               + k1 * num_mol_interp_pts_per_node * num_mol_interp_pts_per_node
                               + k2 * num_mol_interp_pts_per_node + k3;
                               
                double dist_x = xx - mol_clusters_x_ptr[node_interp_pts_start + k1];
                double dist_y = yy - mol_clusters_y_ptr[node_interp_pts_start + k2];
                double dist_z = zz - mol_clusters_z_ptr[node_interp_pts_start + k3];
                
                double numerator = 1.;

                // If exact_idx == -1, then no issues.
                // If exact_idx != -1, then we want to zero out terms EXCEPT when exactInd=k1.
                if (exact_idx_x == -1) {
                    numerator *= weights_ptr[k1] / dist_x;
                } else {
                    if (exact_idx_x != k1) numerator *= 0.;
                }

                if (exact_idx_y == -1) {
                    numerator *= weights_ptr[k2] / dist_y;
                } else {
                    if (exact_idx_y != k2) numerator *= 0.;
                }

                if (exact_idx_z == -1) {
                    numerator *= weights_ptr[k3] / dist_z;
                } else {
                    if (exact_idx_z != k3) numerator *= 0.;
                }

                pot_temp += numerator * denominator * mol_clusters_p_ptr[kk];
            }
            }
            }

#ifdef OPENACC_ENABLED
            #pragma acc atomic update
#endif
            coul_eng_ptr[0] += pot_temp * qq;
        }
    } //end loop over nodes
#ifdef OPENACC_ENABLED
    #pragma acc exit data delete(weights_ptr[0:weights_num])
#endif

//    timers_.downward_pass.stop();
*/
}



void CoulombicEnergyCompute::copyin_clusters_to_device() const
{
//    timers_.copyin_clusters_to_device.start();

#ifdef OPENACC_ENABLED
    const double* q_ptr = mol_interp_charge_.data();
    std::size_t q_num   = mol_interp_charge_.size();
    
    const double* p_ptr = mol_interp_potential_.data();
    std::size_t p_num   = mol_interp_potential_.size();

    const double* coul_eng_ptr = coul_eng_vec_.data();
    std::size_t coul_eng_num   = coul_eng_vec_.size();
    
    #pragma acc enter data copyin(q_ptr[0:q_num], p_ptr[0:p_num])
    #pragma acc enter data copyin(coul_eng_ptr[0:coul_eng_num])
#endif

//    timers_.copyin_clusters_to_device.stop();
}


void CoulombicEnergyCompute::delete_clusters_from_device() const
{
//    timers_.delete_clusters_from_device.start();

#ifdef OPENACC_ENABLED
    const double* q_ptr = mol_interp_charge_.data();
    std::size_t q_num   = mol_interp_charge_.size();
    
    const double* p_ptr = mol_interp_potential_.data();
    std::size_t p_num   = mol_interp_potential_.size();
    
    const double* coul_eng_ptr = coul_eng_vec_.data();
    std::size_t coul_eng_num   = coul_eng_vec_.size();

    #pragma acc exit data delete(q_ptr[0:q_num], p_ptr[0:p_num])
    #pragma acc exit data copyout(coul_eng_ptr[0:coul_eng_num])
#endif

//    timers_.delete_clusters_from_device.stop();
}
