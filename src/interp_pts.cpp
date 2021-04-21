#include <cmath>
#include <cstddef>

#include "interp_pts.h"
#include "tree.h"
#include "constants.h"

InterpolationPoints::InterpolationPoints(const class Tree& tree, int degree)
    : tree_(tree)
{
    //timers_.ctor.start();

    num_interp_pts_per_node_ = degree + 1;
    num_interp_pts_ = tree_.num_nodes() * num_interp_pts_per_node_;

    interp_x_.resize(num_interp_pts_);
    interp_y_.resize(num_interp_pts_);
    interp_z_.resize(num_interp_pts_);

    //timers_.ctor.stop();
}


void InterpolationPoints::compute_all_interp_pts()
{
    //timers_.compute_all_interp_pts.start();

    double* __restrict clusters_x_ptr   = interp_x_.data();
    double* __restrict clusters_y_ptr   = interp_y_.data();
    double* __restrict clusters_z_ptr   = interp_z_.data();
    
    int num_interp_pts_per_node = num_interp_pts_per_node_;
    int degree = num_interp_pts_per_node - 1;
    for (std::size_t node_idx = 0; node_idx < tree_.num_nodes(); ++node_idx) {
    
        std::size_t node_start = node_idx * num_interp_pts_per_node_;
        auto node_bounds = tree_.node_particle_bounds(node_idx);

#ifdef OPENACC_ENABLED
        #pragma acc parallel loop present(clusters_x_ptr, clusters_y_ptr, clusters_z_ptr)
#endif
        for (int i = 0; i < num_interp_pts_per_node; ++i) {
            double tt = std::cos(i * constants::PI / degree);
            clusters_x_ptr[node_start + i] = node_bounds[0] + (tt + 1.) / 2. * (node_bounds[1] - node_bounds[0]);
            clusters_y_ptr[node_start + i] = node_bounds[2] + (tt + 1.) / 2. * (node_bounds[3] - node_bounds[2]);
            clusters_z_ptr[node_start + i] = node_bounds[4] + (tt + 1.) / 2. * (node_bounds[5] - node_bounds[4]);
        }
    }

    //timers_.compute_all_interp_pts.stop();
}




void InterpolationPoints::copyin_to_device() const
{
//    timers_.copyin_to_device.start();

#ifdef OPENACC_ENABLED
    const double* x_ptr = interp_x_.data();
    const double* y_ptr = interp_y_.data();
    const double* z_ptr = interp_z_.data();
    
    std::size_t x_num = interp_x_.size();
    std::size_t y_num = interp_y_.size();
    std::size_t z_num = interp_z_.size();
    
    #pragma acc enter data create(x_ptr[0:x_num], y_ptr[0:y_num], z_ptr[0:z_num])
#endif

//    timers_.copyin_to_device.stop();
}


void InterpolationPoints::delete_from_device() const
{
//    timers_.delete_from_device.start();

#ifdef OPENACC_ENABLED
    const double* x_ptr = interp_x_.data();
    const double* y_ptr = interp_y_.data();
    const double* z_ptr = interp_z_.data();
    
    std::size_t x_num = interp_x_.size();
    std::size_t y_num = interp_y_.size();
    std::size_t z_num = interp_z_.size();
    
    #pragma acc exit data delete(x_ptr[0:x_num], y_ptr[0:y_num], z_ptr[0:z_num])
#endif

//    timers_.delete_from_device.stop();
}
