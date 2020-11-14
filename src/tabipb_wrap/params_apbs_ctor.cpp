#include <cmath>
#include <cstdlib>

#include "../constants.h"
#include "../params.h"

Params::Params(TABIPBInput tabipbIn)
{
    if (tabipbIn.mesh_flag_ == 1) mesh_ = SES;
    if (tabipbIn.mesh_flag_ == 2) mesh_ = SKIN;

    mesh_density_ = tabipbIn.mesh_density_;
    mesh_probe_radius_ = tabipbIn.mesh_probe_radius_;
    
    phys_temp_ = tabipbIn.phys_temp_;
    phys_eps_solute_ = tabipbIn.phys_eps_solute_;
    phys_eps_solvent_ = tabipbIn.phys_eps_solvent_;
    phys_bulk_strength_ = tabipbIn.phys_bulk_strength_;
    
    tree_degree_ = tabipbIn.tree_degree_;
    tree_max_per_leaf_ = tabipbIn.tree_max_per_leaf_;
    tree_theta_ = tabipbIn.tree_theta_;

    if (tabipbIn.precondition_ == 1) precondition_ = true;
    if (tabipbIn.nonpolar_ == 1) nonpolar_ = true;
    
    if (tabipbIn.output_data_ == 1) output_vtk_ = true;
    output_csv_ = false;
    output_timers_ = false;
    
    phys_eps_    = phys_eps_solvent_ / phys_eps_solute_;
    phys_kappa2_ = constants::BULK_COEFF * phys_bulk_strength_ / phys_eps_solvent_ / phys_temp_;
    phys_kappa_  = std::sqrt(phys_kappa2_);
}
