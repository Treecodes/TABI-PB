#include <vector>
#include <string>
#include <fstream>
#include <sstream>
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


double Molecule::compute_coulombic_energy() const
{
    double coulombic_energy = 0.0;
    double epsp = params_.phys_eps_solute_;
    
    const double* __restrict__ molecule_coords_ptr = coords_.data();
    const double* __restrict__ molecule_charge_ptr = charge_.data();

    for (std::size_t i = 0; i < num_atoms_; ++i) {
        double i_pos_x  = molecule_coords_ptr[3*i];
        double i_pos_y  = molecule_coords_ptr[3*i + 1];
        double i_pos_z  = molecule_coords_ptr[3*i + 2];
        double i_charge = molecule_charge_ptr[i];
        
        for (std::size_t j = i+1; j < num_atoms_; ++j) {
            double dist_x = i_pos_x - molecule_coords_ptr[3*j];
            double dist_y = i_pos_y - molecule_coords_ptr[3*j + 1];
            double dist_z = i_pos_z - molecule_coords_ptr[3*j + 2];
            coulombic_energy += i_charge * molecule_charge_ptr[j] / epsp
                              / std::sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
        }
    }

    return coulombic_energy;
}
