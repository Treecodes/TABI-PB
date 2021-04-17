#include "../tabipb_timers.h"
#include "../params.h"
#include "../molecule.h"
#include "../elements.h"
#include "../tree.h"
#include "../interaction_list.h"
#include "../clusters.h"
#include "../boundary_element.h"

#include "TABIPBWrap.h"

struct TABIPBOutput runTABIPBWrapAPBS(struct TABIPBInput tabipbIn, Valist* APBSMolecule)
{
    struct Params params(tabipbIn);
    struct Timers timers;

    class Molecule molecule(APBSMolecule, params, timers.molecule);
    
    molecule.compute_coulombic_energy();
    molecule.build_xyzr_file();
    
    class Elements elements(molecule, params, timers.elements);
    class Tree tree(elements, params, timers.tree);
    
    elements.compute_source_term();
    
    class Clusters clusters(elements, tree, params, timers.clusters);
    
    clusters.compute_all_interp_pts();
    
    class InteractionList interaction_list(tree, params, timers.interaction_list);
    class BoundaryElement boundary_element(elements, clusters, tree, interaction_list, molecule,
                            params, timers.boundary_element);
    
    boundary_element.run_GMRES();
    boundary_element.finalize();
    auto energies = Output(boundary_element, timers);

    return TABIPBOutput{energies[0], energies[1], energies[2]};
}
