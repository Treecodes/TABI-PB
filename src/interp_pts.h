#ifndef H_TABIPB_INTERP_PTS_STRUCT_H
#define H_TABIPB_INTERP_PTS_STRUCT_H

#include <cstddef>

#include "tree.h"

class InterpolationPoints
{
private:
    const class Tree& tree_;

    int num_interp_pts_per_node_;
    std::size_t num_interp_pts_;

    std::vector<double> interp_x_;
    std::vector<double> interp_y_;
    std::vector<double> interp_z_;
    
    
public:
    InterpolationPoints(const class Tree&, int degree);
    ~InterpolationPoints() = default;
    
    std::size_t num_interp_pts_per_node() const { return num_interp_pts_per_node_; };
    
    const std::array<std::size_t, 2> cluster_interp_pts_idxs(std::size_t node_idx) const {
        return std::array<std::size_t, 2> {num_interp_pts_per_node_ *  node_idx,
                                           num_interp_pts_per_node_ * (node_idx + 1)};
    };
    
    const double* interp_x_ptr() const { return interp_x_.data(); };
    const double* interp_y_ptr() const { return interp_y_.data(); };
    const double* interp_z_ptr() const { return interp_z_.data(); };
    
    void compute_all_interp_pts();
    void copyin_to_device() const;
    void delete_from_device() const;
    
};

#endif /* H_TABIPB_INTERP_PTS_STRUCT_H */
