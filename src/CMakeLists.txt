# CXX code for standalone
add_executable(tabipb main.cpp
        params.cpp params.h
        molecule.cpp molecule.h
        particles.cpp particles.h
        tree.cpp tree.h
        clusters.cpp clusters.h
        interaction_list.cpp interaction_list.h
        treecode.cpp gmres.cpp treecode.h
        constants.h)

target_compile_features(tabipb PRIVATE cxx_std_11)
#target_link_libraries(tabipb PRIVATE OpenMP::OpenMP_CXX OpenACC::OpenACC_CXX -acc)
#target_compile_definitions(tabipb PRIVATE OPENACC_ENABLED)
#target_compile_options(tabipb PRIVATE -Minfo=accel)

#Math linking is unnecessary for Windows
if(NOT WIN32)
  target_link_libraries(tabipb PRIVATE m)
endif()

install (TARGETS tabipb DESTINATION bin)
