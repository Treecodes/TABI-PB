# TABIPB
This TABIPB repo serves as a submodule to APBS, and also works as a standalone distribution.

##Build Instructions
To build as an independent executable:
```
mkdir build
cd build
cmake ..
make
```
This creates a tabipb executable located at tabipb/build/bin/tabipb. To run an example, navigate to the examples directory, and run ../build/bin/tabipb.

# Building TABIPB on APBS_Sphinx

This should be similar to the Geoflow build. To invoke TABIPB on Sphinx, you need to use the flag '-DENABLE_TABIPB_SPHINX=ON' when you run cmake. Command:

1. 'cd APBS_SPHINX/plugins/TABIPB/src'
2. 'mkdir build;cd build'
3. 'cmake -DENABLE_TABIPB_SPHINX=ON ..'
4. 'make'
5. 'cp tabipb_sph.so ../..'


