#include <iterator>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstddef>

#include "params.h"
#include "molecule.h"


Molecule::Molecule(struct Params& params, struct Timers_Molecule& timers)
    : params_(params), timers_(timers)
{
    timers_.ctor.start();

    std::string line;
    while (std::getline(params.pqr_file_, line)) {
        
        std::istringstream iss(line);
        std::vector<std::string> tokenized_line{std::istream_iterator<std::string> {iss},
                                                std::istream_iterator<std::string> {} };

        if (tokenized_line[0] == "ATOM") {
            coords_.push_back(std::stod(tokenized_line[5]));
            coords_.push_back(std::stod(tokenized_line[6]));
            coords_.push_back(std::stod(tokenized_line[7]));
            charge_.push_back(std::stod(tokenized_line[8]));
            radius_.push_back(std::stod(tokenized_line[9]));
        }
    }
    
    num_atoms_ = radius_.size();

    timers_.ctor.stop();
}


void Molecule::build_xyzr_file() const
{
    timers_.build_xyzr_file.start();

    std::ofstream xyzr_file("molecule.xyzr");
    
    for (std::size_t i = 0; i < num_atoms_; ++i) {
        xyzr_file << coords_[3*i+0] << " " << coords_[3*i+1] << " "
                  << coords_[3*i+2] << " " << radius_[i] << std::endl;
    }

    xyzr_file.close();

    timers_.build_xyzr_file.stop();
}


void Molecule::compute_coulombic_energy()
{
    timers_.compute_coulombic_energy.start();

    double coulombic_energy = 0.;
    double epsp = params_.phys_eps_solute_;
    std::size_t num_atoms = num_atoms_;
    
    const double* __restrict__ molecule_coords_ptr = coords_.data();
    const double* __restrict__ molecule_charge_ptr = charge_.data();

#ifdef OPENACC_ENABLED
    #pragma acc parallel loop gang present(molecule_coords_ptr, molecule_charge_ptr) \
				   reduction(+:coulombic_energy)
#elif OPENMP_ENABLED
    #pragma omp parallel for reduction(+:coulombic_energy)
#endif
    for (std::size_t i = 0; i < num_atoms; ++i) {
        double i_pos_x  = molecule_coords_ptr[3*i];
        double i_pos_y  = molecule_coords_ptr[3*i + 1];
        double i_pos_z  = molecule_coords_ptr[3*i + 2];
        double i_charge = molecule_charge_ptr[i];
        
#ifdef OPENACC_ENABLED
	#pragma acc loop vector reduction(+:coulombic_energy)
#endif
        for (std::size_t j = i+1; j < num_atoms; ++j) {
            double dist_x = i_pos_x - molecule_coords_ptr[3*j];
            double dist_y = i_pos_y - molecule_coords_ptr[3*j + 1];
            double dist_z = i_pos_z - molecule_coords_ptr[3*j + 2];
            coulombic_energy += i_charge * molecule_charge_ptr[j] / epsp
                              / std::sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
        }
    }

    coulombic_energy_ = coulombic_energy;

    timers_.compute_coulombic_energy.stop();
}


void Molecule::copyin_to_device() const
{
    timers_.copyin_to_device.start();

#ifdef OPENACC_ENABLED
    const double* coords_ptr = coords_.data();
    const double* charge_ptr = charge_.data();

    std::size_t coords_num = coords_.size();
    std::size_t charge_num = charge_.size();

    #pragma acc enter data copyin(coords_ptr[0:coords_num], \
                                  charge_ptr[0:charge_num])
#endif

    timers_.copyin_to_device.stop();
}


void Molecule::delete_from_device() const
{
    timers_.delete_from_device.start();

#ifdef OPENACC_ENABLED
    const double* coords_ptr = coords_.data();
    const double* charge_ptr = charge_.data();

    #pragma acc exit data delete(coords_ptr, charge_ptr)
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
