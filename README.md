# TABIPB
TABI (treecode-accelerated boundary integral) solves the linear Poisson-Boltzmann equation. The solver employs a well-conditioned boundary integral formulation for the electrostatic potential and its normal derivative on the molecular surface, which is triangulated and the integral equations are discretized by centroid collocation. The linear system is solved by GMRES iteration and the matrix-vector product is carried out by a Cartesian terraced which reduces the cost from O(N^2) to O(N\*logN), where N is the number elements.

This TABIPB repo serves as a submodule to APBS, and also works as a standalone distribution. For more information on building APBS, [visit this repository](https://github.com/Electrostatics/apbs-pdb2pqr/tree/master/apbs).


REFERENCE: W.H. Geng and R. Krasny, A treecode-accelerated boundary integral Poisson-Boltzmann solver for continuum electrostatics of solvated biomolecules, J. Comput. Phys. 247, 62-87 (2013)

This material is based upon work supported under NSF Grant DMS-0915057, DMS-1418966, DMS-1418957. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the authors and do not necessarily reflect the views of the NSF.

## Build Instructions
To build as an independent executable:
```
mkdir build; cd build
cmake ..
make
```
This creates a tabipb executable located at tabipb/build/bin/tabipb. To run an example, navigate to the examples directory, and run ../build/bin/tabipb:
```
cd example/
./../build/bin/tabipb
```

## Building TABIPB on APBS_Sphinx

This is similar to the APBS Geoflow build. To invoke TABIPB on Sphinx, you need to use the flag `-DENABLE_TABIPB_SPHINX=ON` when you run cmake. Command:

1. `cd APBS_SPHINX/plugins/TABIPB/src`
2. `mkdir build; cd build`
3. `cmake -DENABLE_TABIPB_SPHINX=ON ..`
4. `make`
5. `cp tabipb_sph.so ../..`


