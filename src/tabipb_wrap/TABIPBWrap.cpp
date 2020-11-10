#include "params.h"
#include "molecule.h"
#include "particles.h"
#include "tree.h"
#include "interaction_list.h"
#include "clusters.h"
#include "treecode.h"

#include "TABIPBWrap.h"

int runTABIPBWrapAPBS(TABIPBInput tabipbIn, Valist* APBSMolecule)
{

    struct Params params(tabipbIn);

    class Molecule molecule(APBSMolecule, params);
    
    molecule.compute_coulombic_energy();
    molecule.build_xyzr_file();
    
    class Particles particles(molecule, params);
    class Tree tree(particles, params);
    
    particles.compute_source_term();
    
    class Clusters clusters(particles, tree, params);
    
    clusters.compute_all_interp_pts();
    
    class InteractionList interaction_list(tree, params);
    class Treecode treecode(particles, clusters,
                            tree, interaction_list, molecule, params);
    
    treecode.run_GMRES();
    treecode.output();

    return 0;
}
