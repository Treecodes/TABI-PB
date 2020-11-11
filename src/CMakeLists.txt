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
if (ENABLE_OPENACC)
    target_link_libraries(tabipb PRIVATE OpenACC::OpenACC_CXX -acc)
    target_compile_definitions(tabipb PRIVATE OPENACC_ENABLED)
    target_compile_options(tabipb PRIVATE -Minfo=accel)
endif ()

if (ENABLE_OPENMP)
    target_link_libraries(tabipb PRIVATE OpenMP::OpenMP_CXX)
endif ()

#Math linking is unnecessary for Windows
if (NOT WIN32)
    target_link_libraries(tabipb PRIVATE m)
endif ()

install (TARGETS tabipb DESTINATION bin)



################################################
###### For APBS 
################################################

#option(ENABLE_TABIPB_APBS "TABIPB APBS" OFF)
if (ENABLE_TABIPB_APBS)
    add_definitions(-DTABIPB_APBS)
    include_directories(
        ../../../src/
        ../../../build/src/
        ../../../externals/fetk/maloc/src/base
        ../../../externals/fetk/maloc/src/vsys
        ../../../externals/fetk/maloc/src/vsh
        ../../../externals/fetk/maloc/src/psh
        )
  
    option(BUILD_SHARED_LIBS "Build shared libraries." OFF)
  
    # FFLAGS depend on the compiler
    get_filename_component (C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)
  
    if (C_COMPILER_NAME MATCHES "gcc.*")
      set (CMAKE_C_FLAGS_RELEASE "-c -O3")
      set (CMAKE_C_FLAGS_DEBUG   "")
    elseif (C_COMPILER_NAME MATCHES "icc.*")
      set (CMAKE_C_FLAGS_RELEASE "-fast -c")
      set (CMAKE_C_FLAGS_DEBUG   "")
    else ()
      message ("CMAKE_C_COMPILER full path: " ${CMAKE_C_COMPILER})
      message ("C compiler: " ${C_COMPILER_NAME})
    endif ()
  
    set(LIBFILES
        params.cpp params.h molecule.cpp molecule.h
        particles.cpp particles.h tree.cpp tree.h
        clusters.cpp clusters.h interaction_list.cpp interaction_list.h
        treecode.cpp gmres.cpp treecode.h constants.h
        tabipb_wrap/TABIPBWrap.cpp tabipb_wrap/TABIPBWrap.h
        tabipb_wrap/TABIPBStruct.h tabipb_wrap/params_apbs_ctor.cpp
        tabipb_wrap/molecule_apbs_ctor.cpp)
  
    set(TABIPB_LIBNAME TABIPBlib)
    add_library(${TABIPB_LIBNAME} ${LIBFILES})
    target_compile_features(${TABIPB_LIBNAME} PRIVATE cxx_std_11)

endif ()
