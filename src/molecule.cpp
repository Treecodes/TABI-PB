#include <numeric>
#include <iterator>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstddef>

#include "molecule.h"


Molecule::Molecule(struct Params& params, struct Timers_Molecule& timers)
    : Particles(params), timers_(timers)
{
    timers_.ctor.start();

    std::string line;
    while (std::getline(params.pqr_file_, line)) {
        
        std::istringstream iss(line);
        std::vector<std::string> tokenized_line{std::istream_iterator<std::string> {iss},
                                                std::istream_iterator<std::string> {} };

        if (tokenized_line[0] == "ATOM") {
            x_.push_back(std::stod(tokenized_line[5]));
            y_.push_back(std::stod(tokenized_line[6]));
            z_.push_back(std::stod(tokenized_line[7]));
            charge_.push_back(std::stod(tokenized_line[8]));
            radius_.push_back(std::stod(tokenized_line[9]));
        }
    }
    
    num_ = radius_.size();
    order_.resize(num_);
    std::iota(order_.begin(), order_.end(), 0);

    timers_.ctor.stop();
}


void Molecule::build_xyzr_file() const
{
    timers_.build_xyzr_file.start();

    std::ofstream xyzr_file("molecule.xyzr");
    
    for (std::size_t i = 0; i < num_; ++i) {
        xyzr_file << x_[i] << " " << y_[i] << " "
                  << z_[i] << " " << radius_[i] << std::endl;
    }

    xyzr_file.close();

    timers_.build_xyzr_file.stop();
}


void Molecule::compute_coulombic_energy()
{
    timers_.compute_coulombic_energy.start();

    double coulombic_energy = 0.;
    double epsp = params_.phys_eps_solute_;
    std::size_t num_atoms = num_;
    
    const double* __restrict molecule_x_ptr = x_.data();
    const double* __restrict molecule_y_ptr = y_.data();
    const double* __restrict molecule_z_ptr = z_.data();
    const double* __restrict molecule_charge_ptr = charge_.data();

#ifdef OPENACC_ENABLED
    #pragma acc parallel loop gang present(molecule_x_ptr, molecule_y_ptr, \
                                           molecule_z_ptr, molecule_charge_ptr) \
                                   reduction(+:coulombic_energy)
#elif OPENMP_ENABLED
    #pragma omp parallel for reduction(+:coulombic_energy)
#endif
    for (std::size_t i = 0; i < num_atoms; ++i) {
        double i_pos_x  = molecule_x_ptr[i];
        double i_pos_y  = molecule_y_ptr[i];
        double i_pos_z  = molecule_z_ptr[i];
        double i_charge = molecule_charge_ptr[i];
        
#ifdef OPENACC_ENABLED
	#pragma acc loop vector reduction(+:coulombic_energy)
#endif
        for (std::size_t j = i+1; j < num_atoms; ++j) {
            double dist_x = i_pos_x - molecule_x_ptr[j];
            double dist_y = i_pos_y - molecule_y_ptr[j];
            double dist_z = i_pos_z - molecule_z_ptr[j];
            coulombic_energy += i_charge * molecule_charge_ptr[j] / epsp
                              / std::sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
        }
    }

    coulombic_energy_ = coulombic_energy;

    timers_.compute_coulombic_energy.stop();
}



//void Molecule::compute_coulombic_energy(class InterpolationPoinsts& interp_pts,
//                                        class Tree& tree,
//                                        class InteractionList& interaction_list)
//{
//    timers_.compute_coulombic_energy.start();
//
//    double coulombic_energy = 0.;
//    double epsp = params_.phys_eps_solute_;
//    std::size_t num_atoms = num_;
//
//    int num_charges_per_node = std::pow(interp_pts.num_interp_pts_per_node(), 3);
//    std::size_t num_charges  = tree.num_nodes() * num_charges_per_node;
//
//    std::vector<double> coulomb_potential(num_atoms);
//    std::vector<double> interp_charge(num_charges);
//    std::vector<double> interp_potential(num_charges);
//
//    coulombic_energy_upward_pass();
//
//    for (std::size_t target_node_idx = 0; target_node_idx < tree_.num_nodes(); ++target_node_idx) {
//
//        for (auto source_node_idx : interaction_list.particle_particle(target_node_idx))
//            coulombic_energy_particle_particle_interact(potential,
//                    tree_.node_particle_idxs(target_node_idx), tree_.node_particle_idxs(source_node_idx));
//
//        for (auto source_node_idx : interaction_list.particle_cluster(target_node_idx))
//            coulombic_energy_particle_cluster_interact(potential,
//                    tree_.node_particle_idxs(target_node_idx), source_node_idx);
//
//        for (auto source_node_idx : interaction_list.cluster_particle(target_node_idx))
//            coulombic_energy_cluster_particle_interact(potential,
//                    target_node_idx, tree_.node_particle_idxs(source_node_idx));
//
//        for (auto source_node_idx : interaction_list.cluster_cluster(target_node_idx))
//            coulombic_energy_cluster_cluster_interact(potential, target_node_idx, source_node_idx);
//    }
//
//    coulombic_engergy_downward_pass(potential);
//

//
//    const double* __restrict molecule_x_ptr = x_.data();
//    const double* __restrict molecule_y_ptr = y_.data();
//    const double* __restrict molecule_z_ptr = z_.data();
//    const double* __restrict molecule_charge_ptr = charge_.data();
//
//#ifdef OPENACC_ENABLED
//    #pragma acc parallel loop gang present(molecule_x_ptr, molecule_y_ptr, \
//                                           molecule_z_ptr, molecule_charge_ptr) \
//                                   reduction(+:coulombic_energy)
//#elif OPENMP_ENABLED
//    #pragma omp parallel for reduction(+:coulombic_energy)
//#endif
//    for (std::size_t i = 0; i < num_atoms; ++i) {
//        double i_pos_x  = molecule_x_ptr[i];
//        double i_pos_y  = molecule_y_ptr[i];
//        double i_pos_z  = molecule_z_ptr[i];
//        double i_charge = molecule_charge_ptr[i];
//
//#ifdef OPENACC_ENABLED
//	#pragma acc loop vector reduction(+:coulombic_energy)
//#endif
//        for (std::size_t j = i+1; j < num_atoms; ++j) {
//            double dist_x = i_pos_x - molecule_x_ptr[j];
//            double dist_y = i_pos_y - molecule_y_ptr[j];
//            double dist_z = i_pos_z - molecule_z_ptr[j];
//            coulombic_energy += i_charge * molecule_charge_ptr[j] / epsp
//                              / std::sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
//        }
//    }
//
//    coulombic_energy_ = coulombic_energy;
//
//    timers_.compute_coulombic_energy.stop();
//}


void Molecule::reorder()
{
    apply_order(order_.begin(), order_.end(), charge_.begin());
    apply_order(order_.begin(), order_.end(), radius_.begin());
}


void Molecule::unorder()
{
    apply_unorder(order_.begin(), order_.end(), x_.begin());
    apply_unorder(order_.begin(), order_.end(), y_.begin());
    apply_unorder(order_.begin(), order_.end(), z_.begin());
    
    apply_unorder(order_.begin(), order_.end(), charge_.begin());
    apply_unorder(order_.begin(), order_.end(), radius_.begin());
}


void Molecule::copyin_to_device() const
{
    timers_.copyin_to_device.start();

#ifdef OPENACC_ENABLED
    const double* x_ptr = x_.data();
    const double* y_ptr = y_.data();
    const double* z_ptr = z_.data();
    const double* charge_ptr = charge_.data();

    std::size_t x_num = x_.size();
    std::size_t y_num = y_.size();
    std::size_t z_num = z_.size();
    std::size_t charge_num = charge_.size();

    #pragma acc enter data copyin(x_ptr[0:x_num], y_ptr[0:y_num], z_ptr[0:z_num] \
                                  charge_ptr[0:charge_num])
#endif

    timers_.copyin_to_device.stop();
}


void Molecule::delete_from_device() const
{
    timers_.delete_from_device.start();

#ifdef OPENACC_ENABLED
    const double* x_ptr = x_.data();
    const double* y_ptr = y_.data();
    const double* z_ptr = z_.data();
    const double* charge_ptr = charge_.data();

    #pragma acc exit data delete(x_ptr, y_ptr, z_ptr, charge_ptr)
#endif

    timers_.delete_from_device.stop();
}


void Timers_Molecule::print() const
{
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout.precision(5);
    std::cout << "|...Molecule function times (s)...." << std::endl;
    std::cout << "|   |...ctor.......................: ";
    std::cout << std::setw(12) << std::right << ctor.elapsed_time() << std::endl;
    std::cout << "|   |...compute_coulombic_energy...: ";
    std::cout << std::setw(12) << std::right << compute_coulombic_energy.elapsed_time() << std::endl;
    std::cout << "|   |...build_xyzr_file............: ";
    std::cout << std::setw(12) << std::right << build_xyzr_file.elapsed_time() << std::endl;
#ifdef OPENACC_ENABLED
    std::cout << "|   |...copyin_to_device...........: ";
    std::cout << std::setw(12) << std::right << copyin_to_device.elapsed_time() << std::endl;
    std::cout << "|   |...delete_from_device.........: ";
    std::cout << std::setw(12) << std::right << delete_from_device.elapsed_time() << std::endl;
#endif
    std::cout << "|" << std::endl;
}


std::string Timers_Molecule::get_durations() const
{
    std::string durations;
    durations.append(std::to_string(ctor                     .elapsed_time())).append(", ");
    durations.append(std::to_string(compute_coulombic_energy .elapsed_time())).append(", ");
    durations.append(std::to_string(build_xyzr_file          .elapsed_time())).append(", ");
    durations.append(std::to_string(copyin_to_device         .elapsed_time())).append(", ");
    durations.append(std::to_string(delete_from_device       .elapsed_time())).append(", ");
    
    return durations;
}


std::string Timers_Molecule::get_headers() const
{
    std::string headers;
    headers.append("Molecule ctor, ");
    headers.append("Molecule compute_coulombic_energy, ");
    headers.append("Molecule build_xyzr_file, ");
    headers.append("Molecule copyin_to_device, ");
    headers.append("Molecule delete_from_device, ");
    
    return headers;
}
