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
    std::vector<double> coeffs_x(num_interp_pts_per_node_);
    std::vector<double> coeffs_y(num_interp_pts_per_node_);
    std::vector<double> coeffs_z(num_interp_pts_per_node_);

    double* weights_ptr = weights.data();
    double* coeffs_x_ptr = coeffs_x.data();
    double* coeffs_y_ptr = coeffs_y.data();
    double* coeffs_z_ptr = coeffs_z.data();

    int weights_num = weights.size();
    int coeffs_x_num = coeffs_x.size();
    int coeffs_y_num = coeffs_y.size();
    int coeffs_z_num = coeffs_z.size();

    int num_interp_pts_per_node = num_interp_pts_per_node_;
    
    for (int i = 0; i < num_interp_pts_per_node_; ++i) {
        weights[i] = ((i % 2 == 0)? 1 : -1);
        if (i == 0 || i == num_interp_pts_per_node_-1) weights[i] = ((i % 2 == 0)? 1 : -1) * 0.5;
    }
    
    for (std::size_t node_idx = 0; node_idx < tree_.num_nodes(); ++node_idx) {
        
        auto particle_idxs = tree_.node_particle_idxs(node_idx);
        std::size_t node_interp_pts_start = node_idx * num_interp_pts_per_node_;
        std::size_t node_charges_start    = node_idx * num_charges_per_node_;
        
#ifdef OPENACC_ENABLED
#pragma acc enter data copyin(weights_ptr[0:weights_num], coeffs_x_ptr[0:coeffs_x_num], \
                              coeffs_y_ptr[0:coeffs_y_num], coeffs_z_ptr[0:coeffs_z_num])
#pragma acc parallel loop present(particles_x_ptr, particles_y_ptr, particles_z_ptr, \
                                  sources_q_ptr, sources_q_dx_ptr, sources_q_dy_ptr, sources_q_dz_ptr, \
                                  clusters_x_ptr, clusters_y_ptr, clusters_z_ptr, \
                                  clusters_q_ptr, clusters_q_dx_ptr, clusters_q_dy_ptr, clusters_q_dz_ptr, \
                                  weights_ptr, coeffs_x_ptr, coeffs_y_ptr, coeffs_z_ptr)
#endif
        for (std::size_t i = particle_idxs[0]; i < particle_idxs[1]; ++i) {
        
            double denominator_x = 0.;
            double denominator_y = 0.;
            double denominator_z = 0.;
            
            int exact_idx_x = -1;
            int exact_idx_y = -1;
            int exact_idx_z = -1;
            
            double xx    = particles_x_ptr[i];
            double yy    = particles_y_ptr[i];
            double zz    = particles_z_ptr[i];
            
            double qq_   = sources_q_ptr   [i];
            double qq_dx = sources_q_dx_ptr[i];
            double qq_dy = sources_q_dy_ptr[i];
            double qq_dz = sources_q_dz_ptr[i];
            
#ifdef OPENACC_ENABLED
            #pragma acc loop reduction(+:denominator_x,denominator_y,denominator_z) \
                             reduction(max:exact_idx_x,exact_idx_y,exact_idx_z)
#endif
            for (int j = 0; j < num_interp_pts_per_node; ++j) {
            
                double dist_x = xx - clusters_x_ptr[node_interp_pts_start + j];
                double dist_y = yy - clusters_y_ptr[node_interp_pts_start + j];
                double dist_z = zz - clusters_z_ptr[node_interp_pts_start + j];

                coeffs_x_ptr[j] = weights_ptr[j] / dist_x;
                coeffs_y_ptr[j] = weights_ptr[j] / dist_y;
                coeffs_z_ptr[j] = weights_ptr[j] / dist_z;
                
                denominator_x += coeffs_x_ptr[j];
                denominator_y += coeffs_y_ptr[j];
                denominator_z += coeffs_z_ptr[j];
                
                if (std::abs(dist_x) < std::numeric_limits<double>::min()) exact_idx_x = j;
                if (std::abs(dist_y) < std::numeric_limits<double>::min()) exact_idx_y = j;
                if (std::abs(dist_z) < std::numeric_limits<double>::min()) exact_idx_z = j;
            }
            
            if (exact_idx_x > -1) {
                for (int id = 0; id < coeffs_x_num; ++id) coeffs_x_ptr[id] = 0.;
                //std::memset(coeffs_x.data(), 0, num_interp_pts_per_node_ * sizeof(double));
                coeffs_x_ptr[exact_idx_x] = 1.;
                denominator_x = 1.;
            }
            
            if (exact_idx_y > -1) {
                for (int id = 0; id < coeffs_y_num; ++id) coeffs_y_ptr[id] = 0.;
                //std::memset(coeffs_y.data(), 0, num_interp_pts_per_node_ * sizeof(double));
                coeffs_y_ptr[exact_idx_y] = 1.;
                denominator_y = 1.;
            }
            
            if (exact_idx_z > -1) {
                for (int id = 0; id < coeffs_z_num; ++id) coeffs_z_ptr[id] = 0.;
                //std::memset(coeffs_z.data(), 0, num_interp_pts_per_node_ * sizeof(double));
                coeffs_z_ptr[exact_idx_z] = 1.;
                denominator_z = 1.;
            }
            
            double denominator = 1. / (denominator_x * denominator_y * denominator_z);
            
#ifdef OPENACC_ENABLED
            #pragma acc loop collapse(3)
#endif
            for (int k1 = 0; k1 < num_interp_pts_per_node; ++k1) {
            for (int k2 = 0; k2 < num_interp_pts_per_node; ++k2) {
            for (int k3 = 0; k3 < num_interp_pts_per_node; ++k3) {
                    
                std::size_t kk = node_charges_start
                               + k1 * num_interp_pts_per_node * num_interp_pts_per_node
                               + k2 * num_interp_pts_per_node + k3;
                               
                double charge_coeff = coeffs_x_ptr[k1] * coeffs_y_ptr[k2]
                                    * coeffs_z_ptr[k3] * denominator;
                                    
                clusters_q_ptr   [kk] += charge_coeff * qq_;
                clusters_q_dx_ptr[kk] += charge_coeff * qq_dx;
                clusters_q_dy_ptr[kk] += charge_coeff * qq_dy;
                clusters_q_dz_ptr[kk] += charge_coeff * qq_dz;
            }
            }
            }
        }
#ifdef OPENACC_ENABLED
#pragma acc exit data delete(weights_ptr[0:weights_num], coeffs_x_ptr[0:coeffs_x_num], \
                             coeffs_y_ptr[0:coeffs_y_num], coeffs_z_ptr[0:coeffs_z_num])
#endif
    }
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
    std::vector<double> coeffs_x(num_interp_pts_per_node_);
    std::vector<double> coeffs_y(num_interp_pts_per_node_);
    std::vector<double> coeffs_z(num_interp_pts_per_node_);

    double* weights_ptr = weights.data();
    double* coeffs_x_ptr = coeffs_x.data();
    double* coeffs_y_ptr = coeffs_y.data();
    double* coeffs_z_ptr = coeffs_z.data();

    int weights_num = weights.size();
    int coeffs_x_num = coeffs_x.size();
    int coeffs_y_num = coeffs_y.size();
    int coeffs_z_num = coeffs_z.size();
    
    std::size_t num_particles = particles_.num();
    int num_interp_pts_per_node = num_interp_pts_per_node_;
    
    for (int i = 0; i < num_interp_pts_per_node_; ++i) {
        weights[i] = ((i % 2 == 0)? 1 : -1);
        if (i == 0 || i == num_interp_pts_per_node_-1) weights[i] = ((i % 2 == 0)? 1 : -1) * 0.5;
    }
    
    for (std::size_t node_idx = 0; node_idx < tree_.num_nodes(); ++node_idx) {
        
        auto particle_idxs = tree_.node_particle_idxs(node_idx);
        std::size_t node_interp_pts_start = node_idx * num_interp_pts_per_node_;
        std::size_t node_potentials_start = node_idx * num_charges_per_node_;

#ifdef OPENACC_ENABLED
#pragma acc enter data copyin(weights_ptr[0:weights_num], coeffs_x_ptr[0:coeffs_x_num], \
                              coeffs_y_ptr[0:coeffs_y_num], coeffs_z_ptr[0:coeffs_z_num])
#pragma acc parallel loop present(particles_x_ptr, particles_y_ptr, particles_z_ptr, \
                                  targets_q_ptr, targets_q_dx_ptr, targets_q_dy_ptr, targets_q_dz_ptr, \
                                  clusters_x_ptr, clusters_y_ptr, clusters_z_ptr, \
                                  clusters_p_ptr, clusters_p_dx_ptr, clusters_p_dy_ptr, clusters_p_dz_ptr, \
                                  potential, weights_ptr, coeffs_x_ptr, coeffs_y_ptr, coeffs_z_ptr)
#endif
        for (std::size_t i = particle_idxs[0]; i < particle_idxs[1]; ++i) {
        
            double denominator_x = 0.;
            double denominator_y = 0.;
            double denominator_z = 0.;
            
            int exact_idx_x = -1;
            int exact_idx_y = -1;
            int exact_idx_z = -1;
            
            double xx = particles_x_ptr[i];
            double yy = particles_y_ptr[i];
            double zz = particles_z_ptr[i];
            
#ifdef OPENACC_ENABLED
            #pragma acc loop reduction(+:denominator_x,denominator_y,denominator_z) \
                             reduction(max:exact_idx_x,exact_idx_y,exact_idx_z)
#endif
            for (int j = 0; j < num_interp_pts_per_node; ++j) {
            
                double dist_x = xx - clusters_x_ptr[node_interp_pts_start + j];
                double dist_y = yy - clusters_y_ptr[node_interp_pts_start + j];
                double dist_z = zz - clusters_z_ptr[node_interp_pts_start + j];

                coeffs_x_ptr[j] = weights_ptr[j] / dist_x;
                coeffs_y_ptr[j] = weights_ptr[j] / dist_y;
                coeffs_z_ptr[j] = weights_ptr[j] / dist_z;
                
                denominator_x += coeffs_x_ptr[j];
                denominator_y += coeffs_y_ptr[j];
                denominator_z += coeffs_z_ptr[j];
                
                if (std::abs(dist_x) < std::numeric_limits<double>::min()) exact_idx_x = j;
                if (std::abs(dist_y) < std::numeric_limits<double>::min()) exact_idx_y = j;
                if (std::abs(dist_z) < std::numeric_limits<double>::min()) exact_idx_z = j;
            }
            
            if (exact_idx_x > -1) {
                for (int id = 0; id < coeffs_x_num; ++id) coeffs_x_ptr[id] = 0.;
                //std::memset(coeffs_x.data(), 0, num_interp_pts_per_node_ * sizeof(double));
                coeffs_x_ptr[exact_idx_x] = 1.;
                denominator_x = 1.;
            }
            
            if (exact_idx_y > -1) {
                for (int id = 0; id < coeffs_y_num; ++id) coeffs_y_ptr[id] = 0.;
                //std::memset(coeffs_y.data(), 0, num_interp_pts_per_node_ * sizeof(double));
                coeffs_y_ptr[exact_idx_y] = 1.;
                denominator_y = 1.;
            }
            
            if (exact_idx_z > -1) {
                for (int id = 0; id < coeffs_z_num; ++id) coeffs_z_ptr[id] = 0.;
                //std::memset(coeffs_z.data(), 0, num_interp_pts_per_node_ * sizeof(double));
                coeffs_z_ptr[exact_idx_z] = 1.;
                denominator_z = 1.;
            }
            
            double denominator = 1. / (denominator_x * denominator_y * denominator_z);
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
                               
                double potential_coeff = coeffs_x_ptr[k1] * coeffs_y_ptr[k2] 
                                       * coeffs_z_ptr[k3] * denominator;
                                    
                pot_comp_   += potential_coeff * clusters_p_ptr   [kk];
                pot_comp_dx += potential_coeff * clusters_p_dx_ptr[kk];
                pot_comp_dy += potential_coeff * clusters_p_dy_ptr[kk];
                pot_comp_dz += potential_coeff * clusters_p_dz_ptr[kk];
            }
            }
            }
            
            potential[i]                 += targets_q_ptr[i]    * pot_comp_;
            potential[i + num_particles] += targets_q_dx_ptr[i] * pot_comp_dx
                                          + targets_q_dy_ptr[i] * pot_comp_dy
                                          + targets_q_dz_ptr[i] * pot_comp_dz;
        }
#ifdef OPENACC_ENABLED
#pragma acc exit data delete(weights_ptr[0:weights_num], coeffs_x_ptr[0:coeffs_x_num], \
                             coeffs_y_ptr[0:coeffs_y_num], coeffs_z_ptr[0:coeffs_z_num])
#endif
    }
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
