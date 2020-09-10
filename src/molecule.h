#ifndef H_TABIPB_MOLECULE_STRUCT_H
#define H_TABIPB_MOLECULE_STRUCT_H

#include <vector>
#include <string>
#include <fstream>
#include <cstddef>

#include "params.h"

class Molecule {

private:
    const struct Params& params_;
    
    std::size_t num_atoms_;
    std::vector<double> coords_;
    std::vector<double> charge_;
    std::vector<double> radius_;

    double coulombic_energy_;

    void compute_coulombic_energy();

public:
    Molecule(struct Params&);
    ~Molecule() = default;
    
    void build_xyzr_file() const;
    
    std::size_t num_atoms() const { return num_atoms_; };
    double coulombic_energy() const { return coulombic_energy_; };
    const double* coords_ptr() const { return coords_.data(); };
    const double* charge_ptr() const { return charge_.data(); };
    const double* radius_ptr() const { return radius_.data(); };
    
    void copyin_to_device() const;
    void delete_from_device() const;
};

#endif /* H_MOLECULE_STRUCT_H */
