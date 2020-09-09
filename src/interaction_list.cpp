#include <iostream>
#include <cmath>
#include <cstddef>

#include "interaction_list.h"

InteractionList::InteractionList(const class Tree& tree, const struct Params& params)
    : tree_(tree), params_(params)
{
    size_check_ = std::pow(params_.tree_degree_ + 1, 3);
    particle_particle_.resize(tree_.num_nodes_);
    particle_cluster_ .resize(tree_.num_nodes_);
    cluster_particle_ .resize(tree_.num_nodes_);
    cluster_cluster_  .resize(tree_.num_nodes_);
    
    //for (auto batch_idx : tree_.leaves_) InteractionList::build_BLTC_lists(batch_idx, 0);
    
    InteractionList::build_BLDTT_lists(0,0);
    
//    for (std::size_t node_idx = 0; node_idx < particle_particle_.size(); ++node_idx) {
//        std::cout << "Node " << node_idx << ", PP: ";
//        for (std::size_t list_idx = 0; list_idx < particle_particle_[node_idx].size(); ++list_idx)
//            std::cout << particle_particle_[node_idx][list_idx] << ", ";
//
//        std::cout << std::endl << "Node " << node_idx << ", PC: ";
//        for (std::size_t list_idx = 0; list_idx < particle_cluster_[node_idx].size(); ++list_idx)
//            std::cout << particle_cluster_[node_idx][list_idx] << ", ";
//
//        std::cout << std::endl << "Node " << node_idx << ", CP: ";
//        for (std::size_t list_idx = 0; list_idx < cluster_particle_[node_idx].size(); ++list_idx)
//            std::cout << cluster_particle_[node_idx][list_idx] << ", ";
//
//        std::cout << std::endl << "Node " << node_idx << ", CC: ";
//        for (std::size_t list_idx = 0; list_idx < cluster_cluster_[node_idx].size(); ++list_idx)
//            std::cout << cluster_cluster_[node_idx][list_idx] << ", ";
//
//        std::cout << std::endl << std::endl;
//    }
}


void InteractionList::build_BLTC_lists(std::size_t batch_idx, std::size_t node_idx)
{
    double dist_x = tree_.node_x_mid_[batch_idx] - tree_.node_x_mid_[node_idx];
    double dist_y = tree_.node_y_mid_[batch_idx] - tree_.node_y_mid_[node_idx];
    double dist_z = tree_.node_z_mid_[batch_idx] - tree_.node_z_mid_[node_idx];
    
    double dist = std::sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
    
    if ((tree_.node_radius_[batch_idx] + tree_.node_radius_[node_idx]) < dist * params_.tree_theta_//) {
       && tree_.node_num_particles_[node_idx] > size_check_) {
       particle_cluster_[batch_idx].push_back(node_idx);
       
    } else if (tree_.node_num_children_[node_idx] == 0) {
        particle_particle_[batch_idx].push_back(node_idx);
    
    } else {
        for (int i = 0; i < tree_.node_num_children_[node_idx]; ++i)
            InteractionList::build_BLTC_lists(batch_idx, tree_.node_children_idx_[8*node_idx + i]);
    }
}


void InteractionList::build_BLDTT_lists(std::size_t target_node_idx, std::size_t source_node_idx)
{
    double dist_x = tree_.node_x_mid_[target_node_idx] - tree_.node_x_mid_[source_node_idx];
    double dist_y = tree_.node_y_mid_[target_node_idx] - tree_.node_y_mid_[source_node_idx];
    double dist_z = tree_.node_z_mid_[target_node_idx] - tree_.node_z_mid_[source_node_idx];
    
    double accept_distance = std::sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z) * params_.tree_theta_;
    double sum_node_radius = tree_.node_radius_[target_node_idx] + tree_.node_radius_[source_node_idx];

    bool target_node_size_check_passed = tree_.node_num_particles_[target_node_idx] > size_check_;
    bool source_node_size_check_passed = tree_.node_num_particles_[source_node_idx] > size_check_;
    
    int target_node_num_children = tree_.node_num_children_[target_node_idx];
    int source_node_num_children = tree_.node_num_children_[source_node_idx];
    
    std::size_t target_node_num_particles = tree_.node_num_particles_[target_node_idx];
    std::size_t source_node_num_particles = tree_.node_num_particles_[source_node_idx];
    
    
    if (sum_node_radius < accept_distance) {
    
        if (!target_node_size_check_passed && !source_node_size_check_passed) {
            particle_particle_[target_node_idx].push_back(source_node_idx);
        
        } else if (!source_node_size_check_passed) {
            cluster_particle_[target_node_idx].push_back(source_node_idx);
            
        } else if (!target_node_size_check_passed) {
            particle_cluster_[target_node_idx].push_back(source_node_idx);
            
        } else {
            cluster_cluster_[target_node_idx].push_back(source_node_idx);
        }
       
    } else {
    
        if (!target_node_num_children && !source_node_num_children) {
            particle_particle_[target_node_idx].push_back(source_node_idx);
    
        } else if (!source_node_num_children) {
            for (int i = 0; i < target_node_num_children; ++i)
                InteractionList::build_BLDTT_lists(tree_.node_children_idx_[8*target_node_idx + i],
                                                   source_node_idx);
    
        } else if (!target_node_num_children) {
            for (int i = 0; i < source_node_num_children; ++i)
                InteractionList::build_BLDTT_lists(target_node_idx,
                                                   tree_.node_children_idx_[8*source_node_idx + i]);
    
        } else if (source_node_num_particles < target_node_num_particles) {
            for (int i = 0; i < target_node_num_children; ++i)
                InteractionList::build_BLDTT_lists(tree_.node_children_idx_[8*target_node_idx + i],
                                                   source_node_idx);
    
        } else {
            for (int i = 0; i < source_node_num_children; ++i)
                InteractionList::build_BLDTT_lists(target_node_idx,
                                                   tree_.node_children_idx_[8*source_node_idx + i]);
        
        }
    }
}
