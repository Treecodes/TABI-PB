#include <vector>
#include <cmath>
#include <cstddef>

#include "../params.h"
#include "../molecule.h"

Molecule::Molecule(Valist* APBSMolecule, 
                   struct Params& params, struct Timers_Molecule& timers)
    : params_(params), timers_(timers)
{
    num_atoms_ = Valist_getNumberAtoms(APBSMolecule);
    
    for (std::size_t i = 0; i < num_atoms_; ++i) {
        Vatom* atom = Valist_getAtom(APBSMolecule, i);
        coords_.push_back(Vatom_getPosition(atom)[0]);
        coords_.push_back(Vatom_getPosition(atom)[1]);
        coords_.push_back(Vatom_getPosition(atom)[2]);
        charge_.push_back(Vatom_getCharge(atom));
        radius_.push_back(Vatom_getRadius(atom));
    }
}
