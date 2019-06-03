# C code for standalone
add_executable(tabipb
    main.c tabipb.c gmres.c
    partition.c run_gmres.c readin.c treecode.c
    print_output.c utilities.c xmalloc.c)

target_compile_definitions(tabipb PUBLIC ${BLA_DEFINITIONS})
target_link_libraries(tabipb ${BLAS_LIBRARIES})
target_include_directories(tabipb
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
    PRIVATE ${BLA_INCLUDE_DIR})

if(ENABLE_MPI_BUILD)
    target_link_libraries(tabipb ${MPI_C_LIBRARIES})
    target_include_directories(tabipb
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        PRIVATE ${MPI_C_INCLUDE_PATH})
endif()

install (TARGETS tabipb DESTINATION bin)
