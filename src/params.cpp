#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

#include "constants.h"
#include "params.h"

Params::Params(char* infile)
{
    std::ifstream paramfile (infile, std::ifstream::in);
    assert(paramfile.good() && "param file is not readable");

    output_dat_ = false;
    output_vtk_ = false;
    output_csv_ = false;
    
    std::string line;
    
    while (std::getline(paramfile, line)) {
        
        std::istringstream iss(line);
        std::vector<std::string> tokenized_line{std::istream_iterator<std::string> {iss},
                                                std::istream_iterator<std::string> {} };
                                        
        std::string param_token = tokenized_line[0];
        std::string param_value = tokenized_line[1];
        
        std::transform(param_token.begin(), param_token.end(), param_token.begin(),
                       [](unsigned char c){ return std::tolower(c); });
                       
        std::transform(param_value.begin(), param_value.end(), param_value.begin(),
                       [](unsigned char c){ return std::tolower(c); });
        
        if (param_token == "mol" || param_token == "pqr") {
            pqr_file_ = std::ifstream(param_value, std::ifstream::in);
            assert(pqr_file_.good() && "pqr file is not readable");
            
        } else if (param_token == "pdie") {
            phys_eps_solute_ = std::stod(param_value);
        
        } else if (param_token == "sdie") {
            phys_eps_solvent_ = std::stod(param_value);
            
        } else if (param_token == "bulk") {
            phys_bulk_strength_ = std::stod(param_value);
        
        } else if (param_token == "temp") {
            phys_temp_ = std::stod(param_value);

        } else if (param_token == "tree_degree") {
            tree_degree_ = std::stoi(param_value);
            assert(tree_degree_ >= 0 && "invalide tree_degree value");

        } else if (param_token == "tree_theta") {
            tree_theta_ = std::stod(param_value);
            assert((tree_theta_ >= 0. && tree_theta_ <= 1.) && "invalide tree_theta value");
        
        } else if (param_token == "tree_max_per_leaf") {
            tree_max_per_leaf_ = std::stoi(param_value);
            assert(tree_max_per_leaf_ >= 0 && "invalid tree_max_per_leaf value");
        
        } else if (param_token == "mesh") {
            auto it = mesh_table_.find(param_value);
            assert(it != mesh_table_.end() && "invalid mesh value");
            mesh_ = it->second;
            
        } else if (param_token == "sdens") {
            mesh_density_ = std::stod(param_value);
            assert(mesh_density_ >= 0 && "invalid density value");
            
        } else if (param_token == "srad") {
            mesh_probe_radius_ = std::stod(param_value);
            assert(mesh_probe_radius_ >= 0 && "invalid probe radius value");
        
        } else if (param_token == "precondition") {
            if (param_value == "true") precondition_ = true;
        
        } else if (param_token == "nonpolar") {
            if (param_value == "true") nonpolar_ = true;
        
        } else if (param_token == "outdata") {
             if (param_value == "dat") output_dat_ = true;
             if (param_value == "vtk") output_vtk_ = true;
             if (param_value == "csv") output_csv_ = true;
        
        } else {
            std::cout << "Skipping undefined token: " << param_token << std::endl;
        }
    }
    
    phys_eps_    = phys_eps_solvent_ / phys_eps_solute_;
    phys_kappa2_ = constants::BULK_COEFF * phys_bulk_strength_ / phys_eps_solvent_ / phys_temp_;
    phys_kappa_  = std::sqrt(phys_kappa2_);
}
