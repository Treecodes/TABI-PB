#include <array>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "tree.h"

Tree::Tree(class Particles& particles, 
           const struct Params& params, struct Timers_Tree& timers)
    : particles_(particles), params_(params), timers_(timers)
{
    timers_.ctor.start();

    num_nodes_     = 0;
    num_leaves_    = 0;
    min_leaf_size_ = std::numeric_limits<std::size_t>::max();
    max_leaf_size_ = std::numeric_limits<std::size_t>::min();
    max_depth_     = 0;

    // tree construction begins with a root on level 0, with no parent
    Tree::construct(0, 0, 0, particles_.num());
    particles_.reorder();
    
    leaves_.resize(num_nodes_);
    std::iota(leaves_.begin(), leaves_.end(), 0);
    auto container_end = std::remove_if(leaves_.begin(), leaves_.end(), [this](std::size_t n)
        {return this->node_num_children_[n] > 0; });
    leaves_.erase(container_end, leaves_.end());

    timers_.ctor.stop();
}


void Tree::construct(std::size_t parent, std::size_t current_level,
                     std::size_t begin,  std::size_t end)
{
    std::size_t node_idx = num_nodes_;
    num_nodes_++;

    if (current_level + 1 > max_depth_) max_depth_ = current_level + 1;
    
    auto bounds = particles_.bounds(begin, end);

    node_particles_begin_.push_back(begin);
    node_particles_end_.push_back(end);
    std::size_t num_particles = end-begin;
    node_num_particles_.push_back(num_particles);
    
    double x_min = bounds[0];
    double x_max = bounds[1];
    double y_min = bounds[2];
    double y_max = bounds[3];
    double z_min = bounds[4];
    double z_max = bounds[5];
    
    double x_len = x_max - x_min;
    double y_len = y_max - y_min;
    double z_len = z_max - z_min;
    
    double radius = std::sqrt(x_len*x_len + y_len*y_len + z_len*z_len) / 2.;
    
    node_x_min_.push_back(x_min);
    node_x_max_.push_back(x_max);
    
    node_y_min_.push_back(y_min);
    node_y_max_.push_back(y_max);
    
    node_z_min_.push_back(z_min);
    node_z_max_.push_back(z_max);
    
    node_x_mid_.push_back((x_min + x_max) / 2.);
    node_y_mid_.push_back((y_min + y_max) / 2.);
    node_z_mid_.push_back((z_min + z_max) / 2.);
    
    node_radius_.push_back(radius);
    
    node_num_children_.push_back(0);
    node_children_idx_.insert(node_children_idx_.end(), {0, 0, 0, 0, 0, 0, 0, 0});
    node_parent_idx_.push_back(parent);
    node_level_.push_back(current_level);
    
    if (num_particles > params_.tree_max_per_leaf_) {
    
        std::array<std::size_t, 16> partitioned_bounds;
        int num_children = particles_.partition_8(begin, end, partitioned_bounds);
        
        int child_ctr = -1;
        for (int i = 0; i < num_children; ++i) {
        
            std::size_t child_begin = partitioned_bounds[2*i + 0];
            std::size_t child_end   = partitioned_bounds[2*i + 1];
            std::size_t child_level = current_level + 1;
                
            if (child_begin < child_end) {
            
                child_ctr++;
                node_num_children_[node_idx]++;
                node_children_idx_[8*node_idx + child_ctr] = num_nodes_;
                Tree::construct(node_idx, child_level, child_begin, child_end);
            }
        }
    } else {
    
        num_leaves_++;
        
        if (num_particles < min_leaf_size_) min_leaf_size_ = num_particles;
        if (num_particles > max_leaf_size_) max_leaf_size_ = num_particles;
    }
}


const std::array<double, 12> Tree::node_particle_bounds(std::size_t node_idx) const
{
    return std::array<double, 12> {node_x_min_[node_idx], node_x_max_[node_idx],
                                   node_y_min_[node_idx], node_y_max_[node_idx],
                                   node_z_min_[node_idx], node_z_max_[node_idx]};
}


const std::array<std::size_t, 2> Tree::node_particle_idxs(std::size_t node_idx) const
{
    return std::array<std::size_t, 2> {node_particles_begin_[node_idx],
                                       node_particles_end_[node_idx]};
}


void Timers_Tree::print() const
{
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout.precision(5);
    std::cout << "|...Tree function times (s)...." << std::endl;
    std::cout << "|   |...ctor.......................: ";
    std::cout << std::setw(12) << std::right << ctor.elapsed_time() << std::endl;
    std::cout << "|" << std::endl;
}


std::string Timers_Tree::get_durations() const
{
    std::string durations;
    durations.append(std::to_string(ctor.elapsed_time())).append(", ");
    
    return durations;
}


std::string Timers_Tree::get_headers() const
{
    std::string headers;
    headers.append("Tree ctor, ");
    
    return headers;
}
