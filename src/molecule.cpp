#include <iterator>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstddef>

#include "params.h"
#include "molecule.h"


Molecule::Molecule(struct Params& params) : params_(params)
{
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
}


void Molecule::build_xyzr_file() const
{
    std::ofstream xyzr_file("molecule.xyzr");
    
    for (std::size_t i = 0; i < num_atoms_; ++i) {
        xyzr_file << coords_[3*i+0] << " " << coords_[3*i+1] << " "
                  << coords_[3*i+2] << " " << radius_[i] << std::endl;
    }

    xyzr_file.close();
}


void Molecule::compute_coulombic_energy()
{
    double coulombic_energy = 0.;
    double epsp = params_.phys_eps_solute_;
    std::size_t num_atoms = num_atoms_;
    
    const double* __restrict__ molecule_coords_ptr = coords_.data();
    const double* __restrict__ molecule_charge_ptr = charge_.data();

#ifdef OPENACC_ENABLED
    #pragma acc parallel loop gang present(molecule_coords_ptr, molecule_charge_ptr) \
				   reduction(+:coulombic_energy)
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
}


void Molecule::copyin_to_device() const
{
#ifdef OPENACC_ENABLED
    #pragma acc enter data copyin(coords_.data()[0:coords_.size()], \
                                  charge_.data()[0:coords_.size()])
#endif
}


void Molecule::delete_from_device() const
{
#ifdef OPENACC_ENABLED
    #pragma acc exit data delete(coords_.data(), charge_.data())
#endif
}
