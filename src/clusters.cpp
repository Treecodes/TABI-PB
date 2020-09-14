#include <array>
#include <limits>
#include <cmath>
#include <cstddef>
#include <cstring>

#include <iostream>

#include "constants.h"
#include "clusters.h"

Clusters::Clusters(const class Particles& particles, const class Tree& tree, const struct Params& params)
    : particles_(particles), tree_(tree), params_(params)
{
    num_interp_pts_per_node_ = params_.tree_degree_ + 1;
    num_charges_per_node_    = std::pow(num_interp_pts_per_node_, 3);
    
    num_interp_pts_ = tree_.num_nodes() * num_interp_pts_per_node_;
    num_charges_    = tree_.num_nodes() * num_charges_per_node_;

    interp_x_.resize(num_interp_pts_);
    interp_y_.resize(num_interp_pts_);
    interp_z_.resize(num_interp_pts_);
    
    interp_charge_.resize(num_charges_);
    interp_charge_dx_.resize(num_charges_);
    interp_charge_dy_.resize(num_charges_);
    interp_charge_dz_.resize(num_charges_);
    
    interp_potential_.resize(num_charges_);
    interp_potential_dx_.resize(num_charges_);
    interp_potential_dy_.resize(num_charges_);
    interp_potential_dz_.resize(num_charges_);
}


void Clusters::compute_all_interp_pts()
{
    double* __restrict__ clusters_x_ptr   = interp_x_.data();
    double* __restrict__ clusters_y_ptr   = interp_y_.data();
    double* __restrict__ clusters_z_ptr   = interp_z_.data();
    
    int degree = params_.tree_degree_;
    int num_interp_pts_per_node = num_interp_pts_per_node_;
    for (std::size_t node_idx = 0; node_idx < tree_.num_nodes(); ++node_idx) {
    
        std::size_t node_start = node_idx * num_interp_pts_per_node_;
        auto node_bounds = tree_.node_particle_bounds(node_idx);

#ifdef OPENACC_ENABLED
#pragma acc parallel loop present(clusters_x_ptr, clusters_y_ptr, clusters_z_ptr)
#endif
        for (std::size_t i = 0; i < num_interp_pts_per_node; ++i) {
            double tt = std::cos(i * constants::PI / degree);
            clusters_x_ptr[node_start + i] = node_bounds[0] + (tt + 1.) / 2. * (node_bounds[1] - node_bounds[0]);
            clusters_y_ptr[node_start + i] = node_bounds[2] + (tt + 1.) / 2. * (node_bounds[3] - node_bounds[2]);
            clusters_z_ptr[node_start + i] = node_bounds[4] + (tt + 1.) / 2. * (node_bounds[5] - node_bounds[4]);
        }
    }
}


void Clusters::upward_pass()
{
    const double* __restrict__ clusters_x_ptr   = interp_x_.data();
    const double* __restrict__ clusters_y_ptr   = interp_y_.data();
    const double* __restrict__ clusters_z_ptr   = interp_z_.data();
    
    double* __restrict__ clusters_q_ptr         = interp_charge_.data();
    double* __restrict__ clusters_q_dx_ptr      = interp_charge_dx_.data();
    double* __restrict__ clusters_q_dy_ptr      = interp_charge_dy_.data();
    double* __restrict__ clusters_q_dz_ptr      = interp_charge_dz_.data();
    
    const double* __restrict__ particles_x_ptr  = particles_.x_ptr();
    const double* __restrict__ particles_y_ptr  = particles_.y_ptr();
    const double* __restrict__ particles_z_ptr  = particles_.z_ptr();
    
    const double* __restrict__ sources_q_ptr    = particles_.source_charge_ptr();
    const double* __restrict__ sources_q_dx_ptr = particles_.source_charge_dx_ptr();
    const double* __restrict__ sources_q_dy_ptr = particles_.source_charge_dy_ptr();
    const double* __restrict__ sources_q_dz_ptr = particles_.source_charge_dz_ptr();
        
    std::vector<double> weights (num_interp_pts_per_node_);
    double* weights_ptr = weights.data();
    int weights_num = weights.size();
    
    for (int i = 0; i < weights_num; ++i) {
        weights[i] = ((i % 2 == 0)? 1 : -1);
        if (i == 0 || i == weights_num-1) weights[i] = ((i % 2 == 0)? 1 : -1) * 0.5;
    }

    int num_interp_pts_per_node = num_interp_pts_per_node_;
    
#ifdef OPENACC_ENABLED
    #pragma acc enter data copyin(weights_ptr[0:weights_num])
#endif
    
    for (std::size_t node_idx = 0; node_idx < tree_.num_nodes(); ++node_idx) {
        
        auto particle_idxs = tree_.node_particle_idxs(node_idx);
        
        std::size_t node_interp_pts_start = node_idx * num_interp_pts_per_node_;
        std::size_t node_charges_start    = node_idx * num_charges_per_node_;
        
        std::size_t particle_start = particle_idxs[0];
        std::size_t particle_end   = particle_idxs[1];
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
#pragma acc parallel present(particles_x_ptr, particles_y_ptr, particles_z_ptr, \
                             sources_q_ptr, sources_q_dx_ptr, sources_q_dy_ptr, sources_q_dz_ptr, \
                             clusters_x_ptr, clusters_y_ptr, clusters_z_ptr, \
                             clusters_q_ptr, clusters_q_dx_ptr, clusters_q_dy_ptr, clusters_q_dz_ptr, \
                             weights_ptr) \
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
        #pragma acc loop
#endif
        for (std::size_t i = 0; i < num_particles; ++i) {
        
            double denominator_x = 0.;
            double denominator_y = 0.;
            double denominator_z = 0.;
            
            double xx    = particles_x_ptr[particle_start + i];
            double yy    = particles_y_ptr[particle_start + i];
            double zz    = particles_z_ptr[particle_start + i];
            
#ifdef OPENACC_ENABLED
            #pragma acc loop reduction(+:denominator_x,denominator_y,denominator_z) \
                             reduction(max:exact_idx_x_ptr[i],exact_idx_y_ptr[i],exact_idx_z_ptr[i])
#endif
            for (int j = 0; j < num_interp_pts_per_node; ++j) {
            
                double dist_x = xx - clusters_x_ptr[node_interp_pts_start + j];
                double dist_y = yy - clusters_y_ptr[node_interp_pts_start + j];
                double dist_z = zz - clusters_z_ptr[node_interp_pts_start + j];
                
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
        #pragma acc loop
#endif
        for (int k1 = 0; k1 < num_interp_pts_per_node; ++k1) {
        for (int k2 = 0; k2 < num_interp_pts_per_node; ++k2) {
        for (int k3 = 0; k3 < num_interp_pts_per_node; ++k3) {
        
            std::size_t kk = node_charges_start
                   + k1 * num_interp_pts_per_node * num_interp_pts_per_node
                   + k2 * num_interp_pts_per_node + k3;
                   
            double cx = clusters_x_ptr[node_interp_pts_start + k1];
            double w1 = weights_ptr[k1];

            double cy = clusters_y_ptr[node_interp_pts_start + k2];
            double w2 = weights_ptr[k2];
            
            double cz = clusters_z_ptr[node_interp_pts_start + k3];
            double w3 = weights_ptr[k3];
            
            double q_temp    = 0.;
            double q_dx_temp = 0.;
            double q_dy_temp = 0.;
            double q_dz_temp = 0.;
            
#ifdef OPENACC_ENABLED
            #pragma acc loop vector(32) reduction(+:q_temp,q_dx_temp,q_dy_temp,q_dz_temp)
#endif
            for (int i = 0; i < num_particles; i++) {  // loop over source points
            
                double dist_x = particles_x_ptr[particle_start + i] - cx;
                double dist_y = particles_y_ptr[particle_start + i] - cy;
                double dist_z = particles_z_ptr[particle_start + i] - cz;
                
                double numerator = 1.;

                // If exact_idx[i] == -1, then no issues.
                // If exact_idx[i] != -1, then we want to zero out terms EXCEPT when exactInd=k1.
                if (exact_idx_x[i] == -1) {
                    numerator *= w1 / dist_x;
                } else {
                    if (exact_idx_x[i] != k1) numerator *= 0.;
                }

                if (exact_idx_y[i] == -1) {
                    numerator *= w2 / dist_y;
                } else {
                    if (exact_idx_y[i] != k2) numerator *= 0.;
                }

                if (exact_idx_z[i] == -1) {
                    numerator *= w3 / dist_z;
                } else {
                    if (exact_idx_z[i] != k3) numerator *= 0.;
                }

                q_temp    += sources_q_ptr   [particle_start + i] * numerator * denominator_ptr[i];
                q_dx_temp += sources_q_dx_ptr[particle_start + i] * numerator * denominator_ptr[i];
                q_dy_temp += sources_q_dy_ptr[particle_start + i] * numerator * denominator_ptr[i];
                q_dz_temp += sources_q_dz_ptr[particle_start + i] * numerator * denominator_ptr[i];
            }
            
            clusters_q_ptr   [kk] += q_temp;
            clusters_q_dx_ptr[kk] += q_dx_temp;
            clusters_q_dy_ptr[kk] += q_dy_temp;
            clusters_q_dz_ptr[kk] += q_dz_temp;
        }
        }
        }
        
        } // end parallel region
    } // end loop over nodes
#ifdef OPENACC_ENABLED
    #pragma acc exit data delete(weights_ptr[0:weights_num])
#endif
}


void Clusters::downward_pass(double* __restrict__ potential)
{
    const double* __restrict__ clusters_x_ptr    = interp_x_.data();
    const double* __restrict__ clusters_y_ptr    = interp_y_.data();
    const double* __restrict__ clusters_z_ptr    = interp_z_.data();
    
    const double* __restrict__ clusters_p_ptr    = interp_potential_.data();
    const double* __restrict__ clusters_p_dx_ptr = interp_potential_dx_.data();
    const double* __restrict__ clusters_p_dy_ptr = interp_potential_dy_.data();
    const double* __restrict__ clusters_p_dz_ptr = interp_potential_dz_.data();
    
    const double* __restrict__ particles_x_ptr   = particles_.x_ptr();
    const double* __restrict__ particles_y_ptr   = particles_.y_ptr();
    const double* __restrict__ particles_z_ptr   = particles_.z_ptr();
    
    const double* __restrict__ targets_q_ptr     = particles_.target_charge_ptr();
    const double* __restrict__ targets_q_dx_ptr  = particles_.target_charge_dx_ptr();
    const double* __restrict__ targets_q_dy_ptr  = particles_.target_charge_dy_ptr();
    const double* __restrict__ targets_q_dz_ptr  = particles_.target_charge_dz_ptr();
    
    std::vector<double> weights (num_interp_pts_per_node_);
    double* weights_ptr = weights.data();
    int weights_num = weights.size();
    
    for (int i = 0; i < weights_num; ++i) {
        weights[i] = ((i % 2 == 0)? 1 : -1);
        if (i == 0 || i == weights_num-1) weights[i] = ((i % 2 == 0)? 1 : -1) * 0.5;
    }
    
    std::size_t potential_offset = particles_.num();
    int num_interp_pts_per_node = num_interp_pts_per_node_;
    
#ifdef OPENACC_ENABLED
#pragma acc enter data copyin(weights_ptr[0:weights_num])
#endif
    
    for (std::size_t node_idx = 0; node_idx < tree_.num_nodes(); ++node_idx) {
        
        auto particle_idxs = tree_.node_particle_idxs(node_idx);
        std::size_t node_interp_pts_start = node_idx * num_interp_pts_per_node_;
        std::size_t node_potentials_start = node_idx * num_charges_per_node_;
        
        std::size_t particle_start = particle_idxs[0];
        std::size_t particle_end   = particle_idxs[1];
        std::size_t num_particles  = particle_idxs[1] - particle_idxs[0];

#ifdef OPENACC_ENABLED
#pragma acc parallel loop present(particles_x_ptr, particles_y_ptr, particles_z_ptr, \
                                  targets_q_ptr, targets_q_dx_ptr, targets_q_dy_ptr, targets_q_dz_ptr, \
                                  clusters_x_ptr, clusters_y_ptr, clusters_z_ptr, \
                                  clusters_p_ptr, clusters_p_dx_ptr, clusters_p_dy_ptr, clusters_p_dz_ptr, \
                                  potential, weights_ptr)
#endif
        for (std::size_t i = 0; i < num_particles; ++i) {
        
            double denominator_x = 0.;
            double denominator_y = 0.;
            double denominator_z = 0.;
            
            int exact_idx_x = -1;
            int exact_idx_y = -1;
            int exact_idx_z = -1;
            
            double xx = particles_x_ptr[particle_start + i];
            double yy = particles_y_ptr[particle_start + i];
            double zz = particles_z_ptr[particle_start + i];
            
#ifdef OPENACC_ENABLED
            #pragma acc loop reduction(+:denominator_x,denominator_y,denominator_z) \
                             reduction(max:exact_idx_x,exact_idx_y,exact_idx_z)
#endif
            for (int j = 0; j < num_interp_pts_per_node; ++j) {
            
                double dist_x = xx - clusters_x_ptr[node_interp_pts_start + j];
                double dist_y = yy - clusters_y_ptr[node_interp_pts_start + j];
                double dist_z = zz - clusters_z_ptr[node_interp_pts_start + j];
                
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

            double pot_comp_   = 0.;
            double pot_comp_dx = 0.;
            double pot_comp_dy = 0.;
            double pot_comp_dz = 0.;
            
#ifdef OPENACC_ENABLED
            #pragma acc loop collapse(3) reduction(+:pot_comp_,  pot_comp_dx, \
                                                     pot_comp_dy,pot_comp_dz)
#endif
            for (int k1 = 0; k1 < num_interp_pts_per_node; ++k1) {
            for (int k2 = 0; k2 < num_interp_pts_per_node; ++k2) {
            for (int k3 = 0; k3 < num_interp_pts_per_node; ++k3) {
                    
                std::size_t kk = node_potentials_start
                               + k1 * num_interp_pts_per_node * num_interp_pts_per_node
                               + k2 * num_interp_pts_per_node + k3;
                               
                double dist_x = xx - clusters_x_ptr[node_interp_pts_start + k1];
                double dist_y = yy - clusters_y_ptr[node_interp_pts_start + k2];
                double dist_z = zz - clusters_z_ptr[node_interp_pts_start + k3];
                
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

                pot_comp_   += numerator * denominator * clusters_p_ptr   [kk];
                pot_comp_dx += numerator * denominator * clusters_p_dx_ptr[kk];
                pot_comp_dy += numerator * denominator * clusters_p_dy_ptr[kk];
                pot_comp_dz += numerator * denominator * clusters_p_dz_ptr[kk];
            }
            }
            }
            
            double pot_temp_1 = targets_q_ptr   [particle_start + i] * pot_comp_;
            double pot_temp_2 = targets_q_dx_ptr[particle_start + i] * pot_comp_dx
                              + targets_q_dy_ptr[particle_start + i] * pot_comp_dy
                              + targets_q_dz_ptr[particle_start + i] * pot_comp_dz;
#ifdef OPENACC_ENABLED
            #pragma acc atomic update
#endif
            potential[particle_start + i]                    += pot_temp_1;
#ifdef OPENACC_ENABLED
            #pragma acc atomic update
#endif
            potential[particle_start + i + potential_offset] += pot_temp_2;
        }
    } //end loop over nodes
#ifdef OPENACC_ENABLED
    #pragma acc exit data delete(weights_ptr[0:weights_num])
#endif
}


void Clusters::clear_charges()
{
#ifdef OPENACC_ENABLED
    std::size_t num_charges = num_charges_;
    double* __restrict__ clusters_q_ptr    = interp_charge_.data();
    double* __restrict__ clusters_q_dx_ptr = interp_charge_dx_.data();
    double* __restrict__ clusters_q_dy_ptr = interp_charge_dy_.data();
    double* __restrict__ clusters_q_dz_ptr = interp_charge_dz_.data();
    
    #pragma acc parallel loop present(clusters_q_ptr, clusters_q_dx_ptr, \
                                      clusters_q_dy_ptr, clusters_q_dz_ptr)
    for (std::size_t i = 0; i < num_charges; ++i) {
        clusters_q_ptr[i] = 0.;
        clusters_q_dx_ptr[i] = 0.;
        clusters_q_dy_ptr[i] = 0.;
        clusters_q_dz_ptr[i] = 0.;
    }
#else
    std::fill(interp_charge_.begin(),    interp_charge_.end(),    0);
    std::fill(interp_charge_dx_.begin(), interp_charge_dx_.end(), 0);
    std::fill(interp_charge_dy_.begin(), interp_charge_dy_.end(), 0);
    std::fill(interp_charge_dz_.begin(), interp_charge_dz_.end(), 0);
#endif
}


void Clusters::clear_potentials()
{
#ifdef OPENACC_ENABLED
    std::size_t num_potentials = num_charges_;
    double* __restrict__ clusters_p_ptr    = interp_potential_.data();
    double* __restrict__ clusters_p_dx_ptr = interp_potential_dx_.data();
    double* __restrict__ clusters_p_dy_ptr = interp_potential_dy_.data();
    double* __restrict__ clusters_p_dz_ptr = interp_potential_dz_.data();
    
    #pragma acc parallel loop present(clusters_p_ptr, clusters_p_dx_ptr, \
                                      clusters_p_dy_ptr, clusters_p_dz_ptr)
    for (std::size_t i = 0; i < num_potentials; ++i) {
        clusters_p_ptr[i] = 0.;
        clusters_p_dx_ptr[i] = 0.;
        clusters_p_dy_ptr[i] = 0.;
        clusters_p_dz_ptr[i] = 0.;
    }
#else
    std::fill(interp_potential_.begin(),    interp_potential_.end(),    0);
    std::fill(interp_potential_dx_.begin(), interp_potential_dx_.end(), 0);
    std::fill(interp_potential_dy_.begin(), interp_potential_dy_.end(), 0);
    std::fill(interp_potential_dz_.begin(), interp_potential_dz_.end(), 0);
#endif
}


void Clusters::copyin_to_device() const
{
#ifdef OPENACC_ENABLED
    const double* x_ptr = interp_x_.data();
    const double* y_ptr = interp_y_.data();
    const double* z_ptr = interp_z_.data();
    
    std::size_t x_num = interp_x_.size();
    std::size_t y_num = interp_y_.size();
    std::size_t z_num = interp_z_.size();
    
    const double* q_ptr    = interp_charge_.data();
    const double* q_dx_ptr = interp_charge_dx_.data();
    const double* q_dy_ptr = interp_charge_dy_.data();
    const double* q_dz_ptr = interp_charge_dz_.data();
    
    std::size_t q_num    = interp_charge_.size();
    std::size_t q_dx_num = interp_charge_dx_.size();
    std::size_t q_dy_num = interp_charge_dy_.size();
    std::size_t q_dz_num = interp_charge_dz_.size();
    
    const double* p_ptr    = interp_potential_.data();
    const double* p_dx_ptr = interp_potential_dx_.data();
    const double* p_dy_ptr = interp_potential_dy_.data();
    const double* p_dz_ptr = interp_potential_dz_.data();
    
    std::size_t p_num    = interp_potential_.size();
    std::size_t p_dx_num = interp_potential_dx_.size();
    std::size_t p_dy_num = interp_potential_dy_.size();
    std::size_t p_dz_num = interp_potential_dz_.size();
    
    #pragma acc enter data create(x_ptr[0:x_num], y_ptr[0:y_num], z_ptr[0:z_num], \
                q_ptr[0:q_num], q_dx_ptr[0:q_dx_num], q_dy_ptr[0:q_dy_num], q_dz_ptr[0:q_dz_num], \
                p_ptr[0:p_num], p_dx_ptr[0:p_dx_num], p_dy_ptr[0:p_dy_num], p_dz_ptr[0:p_dz_num])
#endif
}


void Clusters::delete_from_device() const
{
#ifdef OPENACC_ENABLED
    const double* x_ptr = interp_x_.data();
    const double* y_ptr = interp_y_.data();
    const double* z_ptr = interp_z_.data();
    
    std::size_t x_num = interp_x_.size();
    std::size_t y_num = interp_y_.size();
    std::size_t z_num = interp_z_.size();
    
    const double* q_ptr    = interp_charge_.data();
    const double* q_dx_ptr = interp_charge_dx_.data();
    const double* q_dy_ptr = interp_charge_dy_.data();
    const double* q_dz_ptr = interp_charge_dz_.data();
    
    std::size_t q_num    = interp_charge_.size();
    std::size_t q_dx_num = interp_charge_dx_.size();
    std::size_t q_dy_num = interp_charge_dy_.size();
    std::size_t q_dz_num = interp_charge_dz_.size();
    
    const double* p_ptr    = interp_potential_.data();
    const double* p_dx_ptr = interp_potential_dx_.data();
    const double* p_dy_ptr = interp_potential_dy_.data();
    const double* p_dz_ptr = interp_potential_dz_.data();
    
    std::size_t p_num    = interp_potential_.size();
    std::size_t p_dx_num = interp_potential_dx_.size();
    std::size_t p_dy_num = interp_potential_dy_.size();
    std::size_t p_dz_num = interp_potential_dz_.size();
    
    #pragma acc exit data delete(x_ptr[0:x_num], y_ptr[0:y_num], z_ptr[0:z_num], \
                q_ptr[0:q_num], q_dx_ptr[0:q_dx_num], q_dy_ptr[0:q_dy_num], q_dz_ptr[0:q_dz_num], \
                p_ptr[0:p_num], p_dx_ptr[0:p_dx_num], p_dy_ptr[0:p_dy_num], p_dz_ptr[0:p_dz_num])
#endif
}


const std::array<std::size_t, 2> Clusters::cluster_interp_pts_idxs(std::size_t node_idx) const
{
    return std::array<std::size_t, 2> {num_interp_pts_per_node_ *  node_idx,
                                       num_interp_pts_per_node_ * (node_idx + 1)};
}


const std::array<std::size_t, 2> Clusters::cluster_charges_idxs(std::size_t node_idx) const
{
    return std::array<std::size_t, 2> {num_charges_per_node_ *  node_idx,
                                       num_charges_per_node_ * (node_idx + 1)};
}
