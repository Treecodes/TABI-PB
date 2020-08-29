# C code for standalone
  add_executable(tabipb
    main.cpp tabipb.cpp gmres.c
    run_gmres.cpp readin.cpp treecode.cpp
    print_output.cpp utilities.cpp xmalloc.cpp
    tree/tree_linked_list.cpp tree/partition.cpp)


target_link_libraries(tabipb ${LAPACK_LIBRARIES})
include_directories(SYSTEM ${LAPACK_INCLUDE_PATH})

  target_compile_definitions(tabipb PUBLIC ${BLA_DEFINITIONS})
#  target_link_libraries(tabipb ${BLAS_LIBRARIES})
#  target_include_directories(tabipb
#    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
#    PRIVATE ${BLA_INCLUDE_DIR})

  if(ENABLE_MPI_BUILD)
    target_link_libraries(tabipb ${MPI_C_LIBRARIES})
    include_directories(SYSTEM ${MPI_C_INCLUDE_PATH})
  endif()

#Math linking is unnecessary for Windows
  if(NOT WIN32)
    target_link_libraries(tabipb m)
  endif()

  install (TARGETS tabipb DESTINATION bin)
