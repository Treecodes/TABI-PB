#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "constants.h"
#include "boundary_element.h"


BoundaryElement::BoundaryElement(class Elements& elements, const class InterpolationPoints& interp_pts,
         const class Tree& tree, const class InteractionList& interaction_list,
         const class Molecule& molecule, 
         const struct Params& params, struct Timers_BoundaryElement& timers)
    : elements_(elements), interp_pts_(interp_pts), tree_(tree),
      interaction_list_(interaction_list), molecule_(molecule), 
      params_(params), timers_(timers)
{
    timers_.ctor.start();

    potential_.assign(2 * elements_.num(), 0.);
    
    num_charges_per_node_ = std::pow(interp_pts_.num_interp_pts_per_node(), 3);
    num_charges_          = tree_.num_nodes() * num_charges_per_node_;
    
    interp_charge_.resize(num_charges_);
    interp_charge_dx_.resize(num_charges_);
    interp_charge_dy_.resize(num_charges_);
    interp_charge_dz_.resize(num_charges_);
    
    interp_potential_.resize(num_charges_);
    interp_potential_dx_.resize(num_charges_);
    interp_potential_dy_.resize(num_charges_);
    interp_potential_dz_.resize(num_charges_);

    timers_.ctor.stop();
}
          
void BoundaryElement::run_GMRES()
{
    timers_.run_GMRES.start();

    long int restrt = 10;
    long int length = 2 * elements_.num();
    long int ldw    = length;
    long int ldh    = restrt + 1;
    
    // These values are modified on return
    residual_       = 1e-4;
    num_iter_       = 100;

    std::vector<double> work_vec(ldw * (restrt + 4));
    std::vector<double> h_vec   (ldh * (restrt + 2));
    
    double* work = work_vec.data();
    double* h    = h_vec.data();

    BoundaryElement::copyin_clusters_to_device();
    
    int err_code = BoundaryElement::gmres_(length, elements_.source_term_ptr(), potential_.data(),
                                    restrt, work, ldw, h, ldh, num_iter_, residual_);

    BoundaryElement::delete_clusters_from_device();

    if (err_code) {
        std::cout << "GMRES error code " << err_code << ". Exiting.";
        std::exit(1);
    }
    
    std::cout << "GMRES completed. " << num_iter_ << " iterations, " << residual_ << " residual.";

    timers_.run_GMRES.stop();
}


void BoundaryElement::matrix_vector(double alpha, const double* __restrict potential_old,
                             double beta,        double* __restrict potential_new)
{
    timers_.matrix_vector.start();

    double potential_coeff_1 = 0.5 * (1. +      params_.phys_eps_);
    double potential_coeff_2 = 0.5 * (1. + 1. / params_.phys_eps_);
    
    std::size_t potential_num = potential_.size();
    auto* potential_temp = (double *)std::malloc(potential_num * sizeof(double));
    std::memcpy(potential_temp, potential_new, potential_num * sizeof(double));
    std::memset(potential_new, 0, potential_num * sizeof(double));

#ifdef OPENACC_ENABLED
    #pragma acc enter data copyin(potential_old[0:potential_num], \
                              potential_new[0:potential_num])
#endif

    BoundaryElement::clear_cluster_charges();
    BoundaryElement::clear_cluster_potentials();

    elements_.compute_charges(potential_old);
    BoundaryElement::upward_pass();

#ifdef OPENMP_ENABLED
    #pragma omp parallel for
#endif
    for (std::size_t target_node_idx = 0; target_node_idx < tree_.num_nodes(); ++target_node_idx) {
        
        for (auto source_node_idx : interaction_list_.particle_particle(target_node_idx))
            BoundaryElement::particle_particle_interact(potential_new, potential_old,
                    tree_.node_particle_idxs(target_node_idx), tree_.node_particle_idxs(source_node_idx));
    
        for (auto source_node_idx : interaction_list_.particle_cluster(target_node_idx))
            BoundaryElement::particle_cluster_interact(potential_new, 
                    tree_.node_particle_idxs(target_node_idx), source_node_idx);
        
        for (auto source_node_idx : interaction_list_.cluster_particle(target_node_idx))
            BoundaryElement::cluster_particle_interact(potential_new, 
                    target_node_idx, tree_.node_particle_idxs(source_node_idx));
        
        for (auto source_node_idx : interaction_list_.cluster_cluster(target_node_idx))
            BoundaryElement::cluster_cluster_interact(potential_new, target_node_idx, source_node_idx);
    }

#ifdef OPENACC_ENABLED
    #pragma acc wait
#endif
    
    BoundaryElement::downward_pass(potential_new);

#ifdef OPENACC_ENABLED
    #pragma acc exit data copyout(potential_old[0:potential_num], \
                              potential_new[0:potential_num])
#endif
    
    for (std::size_t i = 0; i < potential_.size() / 2; ++i)
        potential_new[i] = beta * potential_temp[i]
                + alpha * (potential_coeff_1 * potential_old[i] - potential_new[i]);
                                             
    for (std::size_t i = potential_.size() / 2; i < potential_.size(); ++i)
        potential_new[i] =  beta * potential_temp[i]
                + alpha * (potential_coeff_2 * potential_old[i] - potential_new[i]);
                
    std::free(potential_temp);

    timers_.matrix_vector.stop();
}


void BoundaryElement::particle_particle_interact(double* __restrict potential,
                                          const double* __restrict potential_old,
                                          std::array<std::size_t, 2> target_node_element_idxs,
                                          std::array<std::size_t, 2> source_node_element_idxs)
{
    timers_.particle_particle_interact.start();

    std::size_t target_node_element_begin = target_node_element_idxs[0];
    std::size_t target_node_element_end   = target_node_element_idxs[1];

    std::size_t source_node_element_begin = source_node_element_idxs[0];
    std::size_t source_node_element_end   = source_node_element_idxs[1];
    
    double eps    = params_.phys_eps_;
    double kappa  = params_.phys_kappa_;
    double kappa2 = params_.phys_kappa2_;
    
    const double* __restrict elements_x_ptr    = elements_.x_ptr();
    const double* __restrict elements_y_ptr    = elements_.y_ptr();
    const double* __restrict elements_z_ptr    = elements_.z_ptr();
    
    const double* __restrict elements_nx_ptr   = elements_.nx_ptr();
    const double* __restrict elements_ny_ptr   = elements_.ny_ptr();
    const double* __restrict elements_nz_ptr   = elements_.nz_ptr();

    const double* __restrict elements_area_ptr = elements_.area_ptr();
    
    std::size_t num_elements = elements_.num();

#ifdef OPENACC_ENABLED
    int stream_id = std::rand() % 3;
    #pragma acc parallel loop async(stream_id) present(elements_x_ptr,  elements_y_ptr,  elements_z_ptr, \
                                      elements_nx_ptr, elements_ny_ptr, elements_nz_ptr, \
                                      elements_area_ptr, potential, potential_old)
#endif
    for (std::size_t j = target_node_element_begin; j < target_node_element_end; ++j) {
        
        double target_x = elements_x_ptr[j];
        double target_y = elements_y_ptr[j];
        double target_z = elements_z_ptr[j];
        
        double target_nx = elements_nx_ptr[j];
        double target_ny = elements_ny_ptr[j];
        double target_nz = elements_nz_ptr[j];
        
        double pot_temp_1 = 0.;
        double pot_temp_2 = 0.;

#ifdef OPENACC_ENABLED
        #pragma acc loop reduction(+:pot_temp_1,pot_temp_2)
#endif
        for (std::size_t k = source_node_element_begin; k < source_node_element_end; ++k) {
        
            double source_x = elements_x_ptr[k];
            double source_y = elements_y_ptr[k];
            double source_z = elements_z_ptr[k];
            
            double source_nx = elements_nx_ptr[k];
            double source_ny = elements_ny_ptr[k];
            double source_nz = elements_nz_ptr[k];
            double source_area = elements_area_ptr[k];
            
            double potential_old_0 = potential_old[k];
            double potential_old_1 = potential_old[k + num_elements];
            
            double dist_x = source_x - target_x;
            double dist_y = source_y - target_y;
            double dist_z = source_z - target_z;
            double r = std::sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);
            
            if (r > 0) {
                double one_over_r = 1. / r;
                double G0 = constants::ONE_OVER_4PI * one_over_r;
                double kappa_r = kappa * r;
                double exp_kappa_r = std::exp(-kappa_r);
                double Gk = exp_kappa_r * G0;
                
                double source_cos  = (source_nx * dist_x + source_ny * dist_y + source_nz * dist_z) * one_over_r;
                double target_cos = (target_nx * dist_x + target_ny * dist_y + target_nz * dist_z) * one_over_r;
                
                double tp1 = G0 * one_over_r;
                double tp2 = (1. + kappa_r) * exp_kappa_r;

                double dot_tqsq = source_nx * target_nx + source_ny * target_ny + source_nz * target_nz;
                double G3 = (dot_tqsq - 3. * target_cos * source_cos) * one_over_r * tp1;
                double G4 = tp2 * G3 - kappa2 * target_cos * source_cos * Gk;

                double L1 = source_cos  * tp1 * (1. - tp2 * eps);
                double L2 = G0 - Gk;
                double L3 = G4 - G3;
                double L4 = target_cos * tp1 * (1. - tp2 / eps);
                
                pot_temp_1 += (L1 * potential_old_0 + L2 * potential_old_1) * source_area;
                pot_temp_2 += (L3 * potential_old_0 + L4 * potential_old_1) * source_area;
            }
        }
        
#ifdef OPENACC_ENABLED
        #pragma acc atomic update
#elif OPENMP_ENABLED
        #pragma omp atomic update
#endif
        potential[j]                += pot_temp_1;
#ifdef OPENACC_ENABLED
        #pragma acc atomic update
#elif OPENMP_ENABLED
        #pragma omp atomic update
#endif
        potential[j + num_elements] += pot_temp_2;
    }

    timers_.particle_particle_interact.stop();
}


void BoundaryElement::particle_cluster_interact(double* __restrict potential,
                                         std::array<std::size_t, 2> target_node_element_idxs,
                                         std::size_t source_node_idx)
{
    timers_.particle_cluster_interact.start();

    std::size_t num_elements   = elements_.num();
    int num_interp_pts_per_node = interp_pts_.num_interp_pts_per_node();
    int num_charges_per_node    = num_charges_per_node_;

    std::size_t target_node_element_begin      = target_node_element_idxs[0];
    std::size_t target_node_element_end        = target_node_element_idxs[1];

    std::size_t source_cluster_interp_pts_begin = source_node_idx * num_interp_pts_per_node;
    std::size_t source_cluster_charges_begin    = source_node_idx * num_charges_per_node;
    
    double eps    = params_.phys_eps_;
    double kappa  = params_.phys_kappa_;
    
    const double* __restrict elements_x_ptr   = elements_.x_ptr();
    const double* __restrict elements_y_ptr   = elements_.y_ptr();
    const double* __restrict elements_z_ptr   = elements_.z_ptr();
    
    const double* __restrict targets_q_ptr     = elements_.target_charge_ptr();
    const double* __restrict targets_q_dx_ptr  = elements_.target_charge_dx_ptr();
    const double* __restrict targets_q_dy_ptr  = elements_.target_charge_dy_ptr();
    const double* __restrict targets_q_dz_ptr  = elements_.target_charge_dz_ptr();
    
    const double* __restrict clusters_x_ptr    = interp_pts_.interp_x_ptr();
    const double* __restrict clusters_y_ptr    = interp_pts_.interp_y_ptr();
    const double* __restrict clusters_z_ptr    = interp_pts_.interp_z_ptr();

    const double* __restrict clusters_q_ptr    = interp_charge_.data();
    const double* __restrict clusters_q_dx_ptr = interp_charge_dx_.data();
    const double* __restrict clusters_q_dy_ptr = interp_charge_dy_.data();
    const double* __restrict clusters_q_dz_ptr = interp_charge_dz_.data();
    
#ifdef OPENACC_ENABLED
    int stream_id = std::rand() % 3;
    #pragma acc parallel loop async(stream_id) present(elements_x_ptr, elements_y_ptr, elements_z_ptr, \
                    targets_q_ptr, targets_q_dx_ptr, targets_q_dy_ptr, targets_q_dz_ptr, \
                    clusters_x_ptr, clusters_y_ptr, clusters_z_ptr, \
                    clusters_q_ptr, clusters_q_dx_ptr, clusters_q_dy_ptr, clusters_q_dz_ptr, \
                    potential)
#endif
    for (std::size_t j = target_node_element_begin; j < target_node_element_end; ++j) {

        double target_x = elements_x_ptr[j];
        double target_y = elements_y_ptr[j];
        double target_z = elements_z_ptr[j];
        
        double pot_comp_   = 0.;
        double pot_comp_dx = 0.;
        double pot_comp_dy = 0.;
        double pot_comp_dz = 0.;
        
#ifdef OPENACC_ENABLED
        #pragma acc loop collapse(3) reduction(+:pot_comp_,   pot_comp_dx, \
                                                 pot_comp_dy, pot_comp_dz)
#endif
        for (int k1 = 0; k1 < num_interp_pts_per_node; ++k1) {
        for (int k2 = 0; k2 < num_interp_pts_per_node; ++k2) {
        for (int k3 = 0; k3 < num_interp_pts_per_node; ++k3) {
                
            std::size_t kk = source_cluster_charges_begin
                           + k1 * num_interp_pts_per_node * num_interp_pts_per_node
                           + k2 * num_interp_pts_per_node + k3;

            double dx = target_x - clusters_x_ptr[source_cluster_interp_pts_begin + k1];
            double dy = target_y - clusters_y_ptr[source_cluster_interp_pts_begin + k2];
            double dz = target_z - clusters_z_ptr[source_cluster_interp_pts_begin + k3];

            double r2    = dx*dx + dy*dy + dz*dz;
            double r     = std::sqrt(r2);
            double rinv  = 1. / r;
            double r3inv = rinv  * rinv * rinv;
            double r5inv = r3inv * rinv * rinv;

            double expkr   =  std::exp(-kappa * r);
            double d1term  =  r3inv * expkr * (1. + (kappa * r));
            double d1term1 = -r3inv + d1term * eps;
            double d1term2 = -r3inv + d1term / eps;
            double d2term  =  r5inv * (-3. + expkr * (3. + (3. * kappa * r)
                                                   + (kappa * kappa * r2)));
            double d3term  =  r3inv * ( 1. - expkr * (1. + kappa * r));

            pot_comp_    += (rinv * (1. - expkr) * (clusters_q_ptr   [kk])
                                      + d1term1 * (clusters_q_dx_ptr[kk] * dx
                                                 + clusters_q_dy_ptr[kk] * dy
                                                 + clusters_q_dz_ptr[kk] * dz));
                                    
            pot_comp_dx  += (clusters_q_ptr   [kk]  * (d1term2 * dx)
                          - (clusters_q_dx_ptr[kk]  * (dx * dx * d2term + d3term)
                          +  clusters_q_dy_ptr[kk]  * (dx * dy * d2term)
                          +  clusters_q_dz_ptr[kk]  * (dx * dz * d2term)));
                         
            pot_comp_dy  += (clusters_q_ptr   [kk]  *  d1term2 * dy
                          - (clusters_q_dx_ptr[kk]  * (dx * dy * d2term)
                          +  clusters_q_dy_ptr[kk]  * (dy * dy * d2term + d3term)
                          +  clusters_q_dz_ptr[kk]  * (dy * dz * d2term)));
                         
            pot_comp_dz  += (clusters_q_ptr   [kk]  *  d1term2 * dz
                          - (clusters_q_dx_ptr[kk]  * (dx * dz * d2term)
                          +  clusters_q_dy_ptr[kk]  * (dy * dz * d2term)
                          +  clusters_q_dz_ptr[kk]  * (dz * dz * d2term + d3term)));
        }
        }
        }
        
#ifdef OPENACC_ENABLED
        #pragma acc atomic update
#elif OPENMP_ENABLED
        #pragma omp atomic update
#endif
        potential[j]                += targets_q_ptr   [j] * pot_comp_;
#ifdef OPENACC_ENABLED
        #pragma acc atomic update
#elif OPENMP_ENABLED
        #pragma omp atomic update
#endif
        potential[j + num_elements] += targets_q_dx_ptr[j] * pot_comp_dx
                                     + targets_q_dy_ptr[j] * pot_comp_dy
                                     + targets_q_dz_ptr[j] * pot_comp_dz;
    }

    timers_.particle_cluster_interact.stop();
}


void BoundaryElement::cluster_particle_interact(double* __restrict potential,
                                         std::size_t target_node_idx,
                                         std::array<std::size_t, 2> source_node_element_idxs)
{
    timers_.cluster_particle_interact.start();

    int num_interp_pts_per_node = interp_pts_.num_interp_pts_per_node();
    int num_potentials_per_node = num_charges_per_node_;
    
    std::size_t target_cluster_interp_pts_begin = target_node_idx * num_interp_pts_per_node;
    std::size_t target_cluster_potentials_begin = target_node_idx * num_potentials_per_node;

    std::size_t source_node_element_begin       = source_node_element_idxs[0];
    std::size_t source_node_element_end         = source_node_element_idxs[1];
    
    double eps    = params_.phys_eps_;
    double kappa  = params_.phys_kappa_;
    
    const double* __restrict clusters_x_ptr    = interp_pts_.interp_x_ptr();
    const double* __restrict clusters_y_ptr    = interp_pts_.interp_y_ptr();
    const double* __restrict clusters_z_ptr    = interp_pts_.interp_z_ptr();
    
    double* __restrict clusters_p_ptr          = interp_potential_.data();
    double* __restrict clusters_p_dx_ptr       = interp_potential_dx_.data();
    double* __restrict clusters_p_dy_ptr       = interp_potential_dy_.data();
    double* __restrict clusters_p_dz_ptr       = interp_potential_dz_.data();
    
    const double* __restrict elements_x_ptr    = elements_.x_ptr();
    const double* __restrict elements_y_ptr    = elements_.y_ptr();
    const double* __restrict elements_z_ptr    = elements_.z_ptr();
    
    const double* __restrict sources_q_ptr     = elements_.source_charge_ptr();
    const double* __restrict sources_q_dx_ptr  = elements_.source_charge_dx_ptr();
    const double* __restrict sources_q_dy_ptr  = elements_.source_charge_dy_ptr();
    const double* __restrict sources_q_dz_ptr  = elements_.source_charge_dz_ptr();
    
#ifdef OPENACC_ENABLED
    int stream_id = std::rand() % 3;
    #pragma acc parallel loop async(stream_id) collapse(3) present(clusters_x_ptr, clusters_y_ptr, clusters_z_ptr, \
                    clusters_p_ptr, clusters_p_dx_ptr, clusters_p_dy_ptr, clusters_p_dz_ptr, \
                    elements_x_ptr, elements_y_ptr, elements_z_ptr, \
                    sources_q_ptr, sources_q_dx_ptr, sources_q_dy_ptr, sources_q_dz_ptr, \
                    potential)
#endif
    for (int j1 = 0; j1 < num_interp_pts_per_node; ++j1) {
    for (int j2 = 0; j2 < num_interp_pts_per_node; ++j2) {
    for (int j3 = 0; j3 < num_interp_pts_per_node; ++j3) {
    
        std::size_t jj = target_cluster_potentials_begin
                       + j1 * num_interp_pts_per_node * num_interp_pts_per_node
                       + j2 * num_interp_pts_per_node + j3;

        double target_x = clusters_x_ptr[target_cluster_interp_pts_begin + j1];
        double target_y = clusters_y_ptr[target_cluster_interp_pts_begin + j2];
        double target_z = clusters_z_ptr[target_cluster_interp_pts_begin + j3];
        
        double pot_comp_   = 0.;
        double pot_comp_dx = 0.;
        double pot_comp_dy = 0.;
        double pot_comp_dz = 0.;
    
#ifdef OPENACC_ENABLED
        #pragma acc loop reduction(+:pot_comp_,   pot_comp_dx, \
                                     pot_comp_dy, pot_comp_dz)
#endif
        for (std::size_t k = source_node_element_begin; k < source_node_element_end; ++k) {

            double dx = target_x - elements_x_ptr[k];
            double dy = target_y - elements_y_ptr[k];
            double dz = target_z - elements_z_ptr[k];

            double r2    = dx*dx + dy*dy + dz*dz;
            double r     = std::sqrt(r2);
            double rinv  = 1. / r;
            double r3inv = rinv  * rinv * rinv;
            double r5inv = r3inv * rinv * rinv;

            double expkr   =  std::exp(-kappa * r);
            double d1term  =  r3inv * expkr * (1. + (kappa * r));
            double d1term1 = -r3inv + d1term * eps;
            double d1term2 = -r3inv + d1term / eps;
            double d2term  =  r5inv * (-3. + expkr * (3. + (3. * kappa * r)
                                                   + (kappa * kappa * r2)));
            double d3term  =  r3inv * ( 1. - expkr * (1. + kappa * r));

            pot_comp_    += (rinv * (1. - expkr) * (sources_q_ptr   [k])
                                      + d1term1 * (sources_q_dx_ptr[k] * dx
                                                 + sources_q_dy_ptr[k] * dy
                                                 + sources_q_dz_ptr[k] * dz));
                                    
            pot_comp_dx  += (sources_q_ptr   [k]  * (d1term2 * dx)
                          - (sources_q_dx_ptr[k]  * (dx * dx * d2term + d3term)
                          +  sources_q_dy_ptr[k]  * (dx * dy * d2term)
                          +  sources_q_dz_ptr[k]  * (dx * dz * d2term)));
                         
            pot_comp_dy  += (sources_q_ptr   [k]  *  d1term2 * dy
                          - (sources_q_dx_ptr[k]  * (dx * dy * d2term)
                          +  sources_q_dy_ptr[k]  * (dy * dy * d2term + d3term)
                          +  sources_q_dz_ptr[k]  * (dy * dz * d2term)));
                         
            pot_comp_dz  += (sources_q_ptr   [k]  *  d1term2 * dz
                          - (sources_q_dx_ptr[k]  * (dx * dz * d2term)
                          +  sources_q_dy_ptr[k]  * (dy * dz * d2term)
                          +  sources_q_dz_ptr[k]  * (dz * dz * d2term + d3term)));
        }
    
#ifdef OPENACC_ENABLED
        #pragma acc atomic update
#elif OPENMP_ENABLED
        #pragma omp atomic update
#endif
        clusters_p_ptr   [jj] += pot_comp_;
#ifdef OPENACC_ENABLED
        #pragma acc atomic update
#elif OPENMP_ENABLED
        #pragma omp atomic update
#endif
        clusters_p_dx_ptr[jj] += pot_comp_dx;
#ifdef OPENACC_ENABLED
        #pragma acc atomic update
#elif OPENMP_ENABLED
        #pragma omp atomic update
#endif
        clusters_p_dy_ptr[jj] += pot_comp_dy;
#ifdef OPENACC_ENABLED
        #pragma acc atomic update
#elif OPENMP_ENABLED
        #pragma omp atomic update
#endif
        clusters_p_dz_ptr[jj] += pot_comp_dz;
    }
    }
    }

    timers_.cluster_particle_interact.stop();
}


void BoundaryElement::cluster_cluster_interact(double* __restrict potential,
                                        std::size_t target_node_idx,
                                        std::size_t source_node_idx)
{
    timers_.cluster_cluster_interact.start();

    int num_interp_pts_per_node = interp_pts_.num_interp_pts_per_node();
    int num_charges_per_node    = num_charges_per_node_;

    std::size_t target_cluster_interp_pts_begin = target_node_idx * num_interp_pts_per_node;
    std::size_t target_cluster_potentials_begin = target_node_idx * num_charges_per_node;
    
    std::size_t source_cluster_interp_pts_begin = source_node_idx * num_interp_pts_per_node;
    std::size_t source_cluster_charges_begin    = source_node_idx * num_charges_per_node;
    
    double eps    = params_.phys_eps_;
    double kappa  = params_.phys_kappa_;
    
    const double* __restrict clusters_x_ptr    = interp_pts_.interp_x_ptr();
    const double* __restrict clusters_y_ptr    = interp_pts_.interp_y_ptr();
    const double* __restrict clusters_z_ptr    = interp_pts_.interp_z_ptr();

    double* __restrict clusters_p_ptr          = interp_potential_.data();
    double* __restrict clusters_p_dx_ptr       = interp_potential_dx_.data();
    double* __restrict clusters_p_dy_ptr       = interp_potential_dy_.data();
    double* __restrict clusters_p_dz_ptr       = interp_potential_dz_.data();
    
    const double* __restrict clusters_q_ptr    = interp_charge_.data();
    const double* __restrict clusters_q_dx_ptr = interp_charge_dx_.data();
    const double* __restrict clusters_q_dy_ptr = interp_charge_dy_.data();
    const double* __restrict clusters_q_dz_ptr = interp_charge_dz_.data();

#ifdef OPENACC_ENABLED
    int stream_id = std::rand() % 3;
    #pragma acc parallel loop async(stream_id) collapse(3) present(clusters_x_ptr, clusters_y_ptr, clusters_z_ptr, \
                    clusters_p_ptr, clusters_p_dx_ptr, clusters_p_dy_ptr, clusters_p_dz_ptr, \
                    clusters_q_ptr, clusters_q_dx_ptr, clusters_q_dy_ptr, clusters_q_dz_ptr, \
                    potential)
#endif
    for (int j1 = 0; j1 < num_interp_pts_per_node; j1++) {
    for (int j2 = 0; j2 < num_interp_pts_per_node; j2++) {
    for (int j3 = 0; j3 < num_interp_pts_per_node; j3++) {
    
        std::size_t jj = target_cluster_potentials_begin
                       + j1 * num_interp_pts_per_node * num_interp_pts_per_node
                       + j2 * num_interp_pts_per_node + j3;

        double target_x = clusters_x_ptr[target_cluster_interp_pts_begin + j1];
        double target_y = clusters_y_ptr[target_cluster_interp_pts_begin + j2];
        double target_z = clusters_z_ptr[target_cluster_interp_pts_begin + j3];
        
        double pot_comp_   = 0.;
        double pot_comp_dx = 0.;
        double pot_comp_dy = 0.;
        double pot_comp_dz = 0.;
    
#ifdef OPENACC_ENABLED
        #pragma acc loop collapse(3) reduction(+:pot_comp_,   pot_comp_dx, \
                                                 pot_comp_dy, pot_comp_dz)
#endif
        for (int k1 = 0; k1 < num_interp_pts_per_node; k1++) {
        for (int k2 = 0; k2 < num_interp_pts_per_node; k2++) {
        for (int k3 = 0; k3 < num_interp_pts_per_node; k3++) {
            
            std::size_t kk = source_cluster_charges_begin
                           + k1 * num_interp_pts_per_node * num_interp_pts_per_node
                           + k2 * num_interp_pts_per_node + k3;

            double dx = target_x - clusters_x_ptr[source_cluster_interp_pts_begin + k1];
            double dy = target_y - clusters_y_ptr[source_cluster_interp_pts_begin + k2];
            double dz = target_z - clusters_z_ptr[source_cluster_interp_pts_begin + k3];

            double r2    = dx*dx + dy*dy + dz*dz;
            double r     = std::sqrt(r2);
            double rinv  = 1.0 / r;
            double r3inv = rinv  * rinv * rinv;
            double r5inv = r3inv * rinv * rinv;

            double expkr   =  std::exp(-kappa * r);
            double d1term  =  r3inv * expkr * (1. + (kappa * r));
            double d1term1 = -r3inv + d1term * eps;
            double d1term2 = -r3inv + d1term / eps;
            double d2term  =  r5inv * (-3. + expkr * (3. + (3. * kappa * r)
                                                   + (kappa * kappa * r2)));
            double d3term  =  r3inv * ( 1. - expkr * (1. + kappa * r));

            pot_comp_    += (rinv * (1. - expkr) * (clusters_q_ptr   [kk])
                                      + d1term1 * (clusters_q_dx_ptr[kk] * dx
                                                 + clusters_q_dy_ptr[kk] * dy
                                                 + clusters_q_dz_ptr[kk] * dz));
                                    
            pot_comp_dx  += (clusters_q_ptr   [kk]  * (d1term2 * dx)
                          - (clusters_q_dx_ptr[kk]  * (dx * dx * d2term + d3term)
                          +  clusters_q_dy_ptr[kk]  * (dx * dy * d2term)
                          +  clusters_q_dz_ptr[kk]  * (dx * dz * d2term)));
                         
            pot_comp_dy  += (clusters_q_ptr   [kk]  *  d1term2 * dy
                          - (clusters_q_dx_ptr[kk]  * (dx * dy * d2term)
                          +  clusters_q_dy_ptr[kk]  * (dy * dy * d2term + d3term)
                          +  clusters_q_dz_ptr[kk]  * (dy * dz * d2term)));
                         
            pot_comp_dz  += (clusters_q_ptr   [kk]  *  d1term2 * dz
                          - (clusters_q_dx_ptr[kk]  * (dx * dz * d2term)
                          +  clusters_q_dy_ptr[kk]  * (dy * dz * d2term)
                          +  clusters_q_dz_ptr[kk]  * (dz * dz * d2term + d3term)));
        }
        }
        }
    
#ifdef OPENACC_ENABLED
        #pragma acc atomic update
#elif OPENMP_ENABLED
        #pragma omp atomic update
#endif
        clusters_p_ptr   [jj] += pot_comp_;
#ifdef OPENACC_ENABLED
        #pragma acc atomic update
#elif OPENMP_ENABLED
        #pragma omp atomic update
#endif
        clusters_p_dx_ptr[jj] += pot_comp_dx;
#ifdef OPENACC_ENABLED
        #pragma acc atomic update
#elif OPENMP_ENABLED
        #pragma omp atomic update
#endif
        clusters_p_dy_ptr[jj] += pot_comp_dy;
#ifdef OPENACC_ENABLED
        #pragma acc atomic update
#elif OPENMP_ENABLED
        #pragma omp atomic update
#endif
        clusters_p_dz_ptr[jj] += pot_comp_dz;
    }
    }
    }

    timers_.cluster_cluster_interact.stop();
}


void BoundaryElement::upward_pass()
{
    timers_.upward_pass.start();

    const double* __restrict clusters_x_ptr   = interp_pts_.interp_x_ptr();
    const double* __restrict clusters_y_ptr   = interp_pts_.interp_y_ptr();
    const double* __restrict clusters_z_ptr   = interp_pts_.interp_z_ptr();
    
    double* __restrict clusters_q_ptr         = interp_charge_.data();
    double* __restrict clusters_q_dx_ptr      = interp_charge_dx_.data();
    double* __restrict clusters_q_dy_ptr      = interp_charge_dy_.data();
    double* __restrict clusters_q_dz_ptr      = interp_charge_dz_.data();
    
    const double* __restrict elements_x_ptr  = elements_.x_ptr();
    const double* __restrict elements_y_ptr  = elements_.y_ptr();
    const double* __restrict elements_z_ptr  = elements_.z_ptr();
    
    const double* __restrict sources_q_ptr    = elements_.source_charge_ptr();
    const double* __restrict sources_q_dx_ptr = elements_.source_charge_dx_ptr();
    const double* __restrict sources_q_dy_ptr = elements_.source_charge_dy_ptr();
    const double* __restrict sources_q_dz_ptr = elements_.source_charge_dz_ptr();
        
    std::vector<double> weights (interp_pts_.num_interp_pts_per_node());
    double* weights_ptr = weights.data();
    int weights_num = weights.size();
    
    for (int i = 0; i < weights_num; ++i) {
        weights[i] = ((i % 2 == 0)? 1 : -1);
        if (i == 0 || i == weights_num-1) weights[i] = ((i % 2 == 0)? 1 : -1) * 0.5;
    }

    int num_interp_pts_per_node = interp_pts_.num_interp_pts_per_node();
    
#ifdef OPENACC_ENABLED
    #pragma acc enter data copyin(weights_ptr[0:weights_num])
#endif
    
    for (std::size_t node_idx = 0; node_idx < tree_.num_nodes(); ++node_idx) {
        
        auto particle_idxs = tree_.node_particle_idxs(node_idx);
        
        std::size_t node_interp_pts_start = node_idx * interp_pts_.num_interp_pts_per_node();
        std::size_t node_charges_start    = node_idx * num_charges_per_node_;
        
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
#pragma acc kernels present(elements_x_ptr, elements_y_ptr, elements_z_ptr, \
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
        #pragma acc loop independent
#endif
        for (std::size_t i = 0; i < num_particles; ++i) {
        
            double denominator_x = 0.;
            double denominator_y = 0.;
            double denominator_z = 0.;
            
            double xx    = elements_x_ptr[particle_start + i];
            double yy    = elements_y_ptr[particle_start + i];
            double zz    = elements_z_ptr[particle_start + i];

            // because there's a reduction over exact_idx[i], this loop carries a
            // backward dependence and won't actually parallelize
#ifdef OPENACC_ENABLED
            #pragma acc loop reduction(+:denominator_x,denominator_y,denominator_z)
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
        #pragma acc loop collapse(3) independent
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
            #pragma acc loop reduction(+:q_temp,q_dx_temp,q_dy_temp,q_dz_temp)
#endif
            for (std::size_t i = 0; i < num_particles; i++) {  // loop over source points
            
                double dist_x = elements_x_ptr[particle_start + i] - cx;
                double dist_y = elements_y_ptr[particle_start + i] - cy;
                double dist_z = elements_z_ptr[particle_start + i] - cz;
                
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

    timers_.upward_pass.stop();
}


void BoundaryElement::downward_pass(double* __restrict potential)
{
    timers_.downward_pass.start();

    const double* __restrict clusters_x_ptr    = interp_pts_.interp_x_ptr();
    const double* __restrict clusters_y_ptr    = interp_pts_.interp_y_ptr();
    const double* __restrict clusters_z_ptr    = interp_pts_.interp_z_ptr();
    
    const double* __restrict clusters_p_ptr    = interp_potential_.data();
    const double* __restrict clusters_p_dx_ptr = interp_potential_dx_.data();
    const double* __restrict clusters_p_dy_ptr = interp_potential_dy_.data();
    const double* __restrict clusters_p_dz_ptr = interp_potential_dz_.data();
    
    const double* __restrict elements_x_ptr   = elements_.x_ptr();
    const double* __restrict elements_y_ptr   = elements_.y_ptr();
    const double* __restrict elements_z_ptr   = elements_.z_ptr();
    
    const double* __restrict targets_q_ptr     = elements_.target_charge_ptr();
    const double* __restrict targets_q_dx_ptr  = elements_.target_charge_dx_ptr();
    const double* __restrict targets_q_dy_ptr  = elements_.target_charge_dy_ptr();
    const double* __restrict targets_q_dz_ptr  = elements_.target_charge_dz_ptr();
    
    std::vector<double> weights (interp_pts_.num_interp_pts_per_node());
    double* weights_ptr = weights.data();
    int weights_num = weights.size();
    
    for (int i = 0; i < weights_num; ++i) {
        weights[i] = ((i % 2 == 0)? 1 : -1);
        if (i == 0 || i == weights_num-1) weights[i] = ((i % 2 == 0)? 1 : -1) * 0.5;
    }
    
    std::size_t potential_offset = elements_.num();
    int num_interp_pts_per_node = interp_pts_.num_interp_pts_per_node();
    
#ifdef OPENACC_ENABLED
#pragma acc enter data copyin(weights_ptr[0:weights_num])
#endif
    
    for (std::size_t node_idx = 0; node_idx < tree_.num_nodes(); ++node_idx) {
        
        auto particle_idxs = tree_.node_particle_idxs(node_idx);
        std::size_t node_interp_pts_start = node_idx * num_interp_pts_per_node;
        std::size_t node_potentials_start = node_idx * num_charges_per_node_;
        
        std::size_t particle_start = particle_idxs[0];
        std::size_t num_particles  = particle_idxs[1] - particle_idxs[0];

#ifdef OPENACC_ENABLED
#pragma acc parallel loop present(elements_x_ptr, elements_y_ptr, elements_z_ptr, \
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
            
            double xx = elements_x_ptr[particle_start + i];
            double yy = elements_y_ptr[particle_start + i];
            double zz = elements_z_ptr[particle_start + i];
            
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

    timers_.downward_pass.stop();
}


void BoundaryElement::clear_cluster_charges()
{
    timers_.clear_cluster_charges.start();

#ifdef OPENACC_ENABLED
    std::size_t num_charges = num_charges_;
    double* __restrict clusters_q_ptr    = interp_charge_.data();
    double* __restrict clusters_q_dx_ptr = interp_charge_dx_.data();
    double* __restrict clusters_q_dy_ptr = interp_charge_dy_.data();
    double* __restrict clusters_q_dz_ptr = interp_charge_dz_.data();
    
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

    timers_.clear_cluster_charges.stop();
}


void BoundaryElement::clear_cluster_potentials()
{
    timers_.clear_cluster_potentials.start();

#ifdef OPENACC_ENABLED
    std::size_t num_potentials = num_charges_;
    double* __restrict clusters_p_ptr    = interp_potential_.data();
    double* __restrict clusters_p_dx_ptr = interp_potential_dx_.data();
    double* __restrict clusters_p_dy_ptr = interp_potential_dy_.data();
    double* __restrict clusters_p_dz_ptr = interp_potential_dz_.data();
    
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

    timers_.clear_cluster_potentials.stop();
}


void BoundaryElement::copyin_clusters_to_device() const
{
    timers_.copyin_clusters_to_device.start();

#ifdef OPENACC_ENABLED
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
    
    #pragma acc enter data create( \
                q_ptr[0:q_num], q_dx_ptr[0:q_dx_num], q_dy_ptr[0:q_dy_num], q_dz_ptr[0:q_dz_num], \
                p_ptr[0:p_num], p_dx_ptr[0:p_dx_num], p_dy_ptr[0:p_dy_num], p_dz_ptr[0:p_dz_num])
#endif

    timers_.copyin_clusters_to_device.stop();
}


void BoundaryElement::delete_clusters_from_device() const
{
    timers_.delete_clusters_from_device.start();

#ifdef OPENACC_ENABLED
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
    
    #pragma acc exit data delete( \
                q_ptr[0:q_num], q_dx_ptr[0:q_dx_num], q_dy_ptr[0:q_dy_num], q_dz_ptr[0:q_dz_num], \
                p_ptr[0:p_num], p_dx_ptr[0:p_dx_num], p_dy_ptr[0:p_dy_num], p_dz_ptr[0:p_dz_num])
#endif

    timers_.delete_clusters_from_device.stop();
}


void BoundaryElement::finalize()
{
    timers_.finalize.start();

    solvation_energy_ = constants::UNITS_PARA  * elements_.compute_solvation_energy(potential_);
    coulombic_energy_ = constants::UNITS_COEFF * molecule_.coulombic_energy();
    free_energy_      = solvation_energy_ + coulombic_energy_;
    
    elements_.unorder(potential_);

    constexpr double pot_scaling = constants::UNITS_COEFF * constants::PI * 4.;
    std::transform(std::begin(potential_), std::end(potential_),
                   std::begin(potential_), [=](double x){ return x * pot_scaling; });
                   
    auto pot_min_max = std::minmax_element(
        potential_.begin(), potential_.begin() + potential_.size() / 2);
                                           
    auto pot_normal_min_max = std::minmax_element(
        potential_.begin() + potential_.size() / 2, potential_.end());
        
    pot_min_ = *pot_min_max.first;
    pot_max_ = *pot_min_max.second;
    
    pot_normal_min_ = *pot_normal_min_max.first;
    pot_normal_max_ = *pot_normal_min_max.second;
    
    timers_.finalize.stop();
}


void Timers_BoundaryElement::print() const
{
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout.precision(5);
    std::cout << "|...BoundaryElement function times (s)...." << std::endl;
    std::cout << "|   |...ctor.......................: ";
    std::cout << std::setw(12) << std::right << ctor                       .elapsed_time() << std::endl;
    std::cout << "|   |...run_GMRES..................: ";
    std::cout << std::setw(12) << std::right << run_GMRES                  .elapsed_time() << std::endl;
    std::cout << "|       |...matrix_vector..........: ";
    std::cout << std::setw(12) << std::right << matrix_vector              .elapsed_time() << std::endl;
    std::cout << "|           |...PP interact........: ";
    std::cout << std::setw(12) << std::right << particle_particle_interact .elapsed_time() << std::endl;
    std::cout << "|           |...PC interact........: ";
    std::cout << std::setw(12) << std::right << particle_cluster_interact  .elapsed_time() << std::endl;
    std::cout << "|           |...CP interact........: ";
    std::cout << std::setw(12) << std::right << cluster_particle_interact  .elapsed_time() << std::endl;
    std::cout << "|           |...CC interact........: ";
    std::cout << std::setw(12) << std::right << cluster_cluster_interact   .elapsed_time() << std::endl;
    std::cout << "|       |...precondition...........: ";
    std::cout << std::setw(12) << std::right << precondition               .elapsed_time() << std::endl;
    std::cout << "|   |...finalize...................: ";
    std::cout << std::setw(12) << std::right << finalize                   .elapsed_time() << std::endl;
    std::cout << "|" << std::endl;
}


std::string Timers_BoundaryElement::get_durations() const
{
    std::string durations;
    durations.append(std::to_string(ctor                       .elapsed_time())).append(", ");
    durations.append(std::to_string(run_GMRES                  .elapsed_time())).append(", ");
    durations.append(std::to_string(matrix_vector              .elapsed_time())).append(", ");
    durations.append(std::to_string(particle_particle_interact .elapsed_time())).append(", ");
    durations.append(std::to_string(particle_cluster_interact  .elapsed_time())).append(", ");
    durations.append(std::to_string(cluster_particle_interact  .elapsed_time())).append(", ");
    durations.append(std::to_string(cluster_cluster_interact   .elapsed_time())).append(", ");
    durations.append(std::to_string(precondition               .elapsed_time())).append(", ");
    durations.append(std::to_string(finalize                   .elapsed_time())).append(", ");
    
    return durations;
}


std::string Timers_BoundaryElement::get_headers() const
{
    std::string headers;
    headers.append("BoundaryElement ctor, ");
    headers.append("BoundaryElement run_GMRES, ");
    headers.append("BoundaryElement matrix_vector, ");
    headers.append("BoundaryElement particle_particle_interact, ");
    headers.append("BoundaryElement particle_cluster_interact, ");
    headers.append("BoundaryElement cluster_particle_interact, ");
    headers.append("BoundaryElement cluster_cluster_interact, ");
    headers.append("BoundaryElement precondition, ");
    headers.append("BoundaryElement finalize, ");
    
    return headers;
}
