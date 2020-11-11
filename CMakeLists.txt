cmake_minimum_required(VERSION 3.12)
project(TABI-PB VERSION 3.0 LANGUAGES CXX C)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Set a default build type if none was specified
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")

  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)

  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

option(ENABLE_OPENACC "OpenACC" OFF)
if(ENABLE_OPENACC)
    find_package(OpenACC REQUIRED)
    set(OPENACC_ACCEL_TARGET tesla)
    add_definitions(-DOPENACC_ENABLED)
endif()

option(ENABLE_OPENMP "OpenMP" OFF)
if(ENABLE_OPENMP)
    find_package(OpenMP)
    add_definitions(-DOPENMP_ENABLED)
endif()

get_filename_component(CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG   "-Wall")

add_subdirectory(src)
