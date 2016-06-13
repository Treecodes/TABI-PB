cmake_minimum_required (VERSION 2.8)
project(tabipb)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(TABIPB_ROOT ${PROJECT_SOURCE_DIR})
set(TABIPB_BUILD ${PROJECT_BINARY_DIR})
set(EXECUTABLE_OUTPUT_PATH ${TABIPB_BUILD}/bin)
set(EXECUTABLE_INSTALL_PATH bin)


################################################
##### Make a directory for the executable
################################################
file(MAKE_DIRECTORY build/bin)

include_directories("/src")

# C code for standalone
add_executable(tabipb
  src/main.c src/gmres.c src/daxpy.c src/dcopy.c src/ddot.c
  src/dgemv.c src/dnrm2.c src/drot.c src/drotg.c src/dscal.c src/d_sign.c
  src/dtrsv.c src/gl_functions.c src/readin.c src/pp_timer.c src/treecode.c
  )


################################################
###### APBS components
################################################

option(ENABLE_TABIPB_APBS "TABIPB APBS" OFF)
if(ENABLE_TABIPB_APBS)

  option(BUILD_SHARED_LIBS "Build shared libraries." OFF)

  # make sure that the default is a RELEASE
  if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE)
  endif (NOT CMAKE_BUILD_TYPE)

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

  set (LIBFILES
    src/tabipb.c src/gmres.c src/daxpy.c src/dcopy.c src/ddot.c
    src/dgemv.c src/dnrm2.c src/drot.c src/drotg.c src/dscal.c src/d_sign.c
    src/dtrsv.c src/gl_functions.c src/readin.c src/pp_timer.c src/treecode.c)

  message(***** ${CMAKE_INSTALL_PREFIX} *****)

  if(NOT DEFINED ${CMAKE_INSTALL_PREFIX})
    set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR})
  endif()

  if(NOT DEFINED CMAKE_INSTALL_LIBDIR)
    set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
  endif()

  set(TABIPB_LIBNAME TABIPBlib)
  add_library(${TABIPB_LIBNAME} ${LIBFILES})

  if(BUILD_SHARED_LIBS)
    add_library(${TABIPB_LIBNAME} ${LIBFILES})
    install(TARGETS ${TABIPB_LIBNAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib${TABIPB_LIBNAME}${CMAKE_SHARED_LIBRARY_SUFFIX}
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
  endif()

endif()

################################################
###### Install
################################################

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
install(TARGETS tabipb RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bin)
