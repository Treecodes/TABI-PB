# CXX code for standalone
add_executable(tabipb main.cpp
        params.cpp params.h
        molecule.cpp molecule.h
        particles.cpp particles.h
        tree.cpp tree.h
        clusters.cpp clusters.h
        interaction_list.cpp interaction_list.h
        treecode.cpp output.cpp 
        gmres.cpp treecode.h
        constants.h)

#target_compile_definitions(tabipb PUBLIC ${BLA_DEFINITIONS})
#target_link_libraries(tabipb ${BLAS_LIBRARIES})
#target_include_directories(tabipb
#    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
#    PRIVATE ${BLA_INCLUDE_DIR})

target_link_libraries(tabipb ${LAPACK_LIBRARIES})
include_directories(SYSTEM ${LAPACK_INCLUDE_PATH})

target_compile_definitions(tabipb PUBLIC ${BLA_DEFINITIONS})
target_compile_features(tabipb PUBLIC cxx_std_11)
#  target_link_libraries(tabipb ${BLAS_LIBRARIES})
#  target_include_directories(tabipb
#    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
#    PRIVATE ${BLA_INCLUDE_DIR})

#Math linking is unnecessary for Windows
if(NOT WIN32)
  target_link_libraries(tabipb m)
endif()

install (TARGETS tabipb DESTINATION bin)
