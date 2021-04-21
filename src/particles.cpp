#include <cmath>

#include "partition.h"
#include "particles.h"


const std::array<double, 6> Particles::bounds(std::size_t begin, std::size_t end) const
{
    auto x_min_max = std::minmax_element(x_.begin() + begin, x_.begin() + end);
    auto y_min_max = std::minmax_element(y_.begin() + begin, y_.begin() + end);
    auto z_min_max = std::minmax_element(z_.begin() + begin, z_.begin() + end);
    
    return std::array<double, 6> {*x_min_max.first, *x_min_max.second,
                                  *y_min_max.first, *y_min_max.second,
                                  *z_min_max.first, *z_min_max.second};
}


int Particles::partition_8(std::size_t begin, std::size_t end,
                           std::array<std::size_t, 16>& partitioned_bounds)
{
    int num_children = 1;
    
    partitioned_bounds[0] = begin;
    partitioned_bounds[1] = end;
    
    auto bounds = Particles::bounds(begin, end);
    
    double x_len = bounds[1] - bounds[0];
    double y_len = bounds[3] - bounds[2];
    double z_len = bounds[5] - bounds[4];
    
    double x_mid = (bounds[1] + bounds[0]) / 2.;
    double y_mid = (bounds[3] + bounds[2]) / 2.;
    double z_mid = (bounds[5] + bounds[4]) / 2.;
    
    double max_len = x_len;
    if (max_len < y_len) max_len = y_len;
    if (max_len < z_len) max_len = z_len;

    double critical_len = max_len / std::sqrt(2.);
    
    bool divide_x = false;
    bool divide_y = false;
    bool divide_z = false;
    
    if (x_len > critical_len) divide_x = true;
    if (y_len > critical_len) divide_y = true;
    if (z_len > critical_len) divide_z = true;

    if (divide_x) {

//  This, unfortunately, does not quite work, but it should be something like this to reorder them in an STL way
//
//        std::vector<size_t> reorder_vec(ind[0][1] - ind[0][0] + 1);
//        std::iota(reorder_vec.begin(), reorder_vec.end(), ind[0][0]);
//
//        auto pivot = std::partition(reorder_vec.begin(), reorder_vec.end(),
//                                    [&x, &x_mid, &ind](size_t elem){ return x[elem] < x_mid; });
//
//        std::transform(reorder_vec.begin(),reorder_vec.end(),reorder_vec.begin(),
//                       [&ind](size_t i){ return i-ind[0][0]; });
//
//        reorder_inplace_destructive(reorder_vec.begin(), reorder_vec.end(),
//                        orderarr.begin()+ind[0][0], x.begin()+ind[0][0], y.begin()+ind[0][0], z.begin()+ind[0][0]);
//
//        ind[1][0] = *pivot + ind[0][0];
//        ind[1][1] = ind[0][1];
//        ind[0][1] = *pivot + ind[0][0] - 1;

        std::size_t node_begin = partitioned_bounds[0];
        std::size_t node_end   = partitioned_bounds[1];

        std::size_t pivot_idx = partition<double>(x_.data(), y_.data(), z_.data(), order_.data(),
                                                  node_begin, node_end, x_mid);
        
        partitioned_bounds[2] = pivot_idx;
        partitioned_bounds[3] = partitioned_bounds[1];
        partitioned_bounds[1] = pivot_idx;

        num_children *= 2;
    }

    if (divide_y) {

        for (int i = 0; i < num_children; ++i) {
        
            std::size_t node_begin = partitioned_bounds[2*i + 0];
            std::size_t node_end   = partitioned_bounds[2*i + 1];
        
            std::size_t pivot_idx = partition<double>(y_.data(), x_.data(), z_.data(), order_.data(),
                                                      node_begin, node_end, y_mid);

            partitioned_bounds[2 * (num_children + i) + 0] = pivot_idx;
            partitioned_bounds[2 * (num_children + i) + 1] = partitioned_bounds[2*i + 1];
            partitioned_bounds[2*i + 1] = pivot_idx;
        }

        num_children *= 2;
    }

    if (divide_z) {

        for (int i = 0; i < num_children; ++i) {
        
            std::size_t node_begin = partitioned_bounds[2*i + 0];
            std::size_t node_end   = partitioned_bounds[2*i + 1];
        
            std::size_t pivot_idx = partition<double>(z_.data(), x_.data(), y_.data(), order_.data(),
                                                      node_begin, node_end, z_mid);

            partitioned_bounds[2 * (num_children + i) + 0] = pivot_idx;
            partitioned_bounds[2 * (num_children + i) + 1] = partitioned_bounds[2*i + 1];
            partitioned_bounds[2*i + 1] = pivot_idx;
        }

        num_children *= 2;
    }

    return num_children;

}
