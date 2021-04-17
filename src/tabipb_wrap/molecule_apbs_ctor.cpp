#include <vector>
#include <cmath>
#include <cstddef>

#include "../params.h"
#include "../molecule.h"

Molecule::Molecule(Valist* APBSMolecule, 
                   struct Params& params, struct Timers_Molecule& timers)
    : params_(params), timers_(timers)
{
    num_ = Valist_getNumberAtoms(APBSMolecule);
    
    for (std::size_t i = 0; i < num_; ++i) {
        Vatom* atom = Valist_getAtom(APBSMolecule, i);
        x_.push_back(Vatom_getPosition(atom)[0]);
        y_.push_back(Vatom_getPosition(atom)[1]);
        z_.push_back(Vatom_getPosition(atom)[2]);
        charge_.push_back(Vatom_getCharge(atom));
        radius_.push_back(Vatom_getRadius(atom));
    }
}
