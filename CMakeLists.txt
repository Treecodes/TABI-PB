cmake_minimum_required(VERSION 3.1...3.14)
project(TABI-PB VERSION 3.0 LANGUAGES CXX C)

#set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Set a default build type if none was specified
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release")
endif()

set(TABIPB_ROOT ${PROJECT_SOURCE_DIR})
set(TABIPB_BUILD ${PROJECT_BINARY_DIR})
set(EXECUTABLE_OUTPUT_PATH ${TABIPB_BUILD}/bin)
set(EXECUTABLE_INSTALL_PATH bin)

################################################
##### Make a directory for the executable
################################################
file(MAKE_DIRECTORY build/bin)

option(ENABLE_OPENACC "OpenACC" OFF)
if(ENABLE_OPENACC)
    add_definitions(-DOPENACC_ENABLED)
    set(OPENACC_ACCEL_TARGET tesla)
    find_package(OpenACC REQUIRED)
endif()
find_package(OpenMP)

get_filename_component(CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG   "-Wall")
message("CMAKE_CXX_COMPILER full path: ${CMAKE_CXX_COMPILER}")
message("CXX compiler: ${C_COMPILER_NAME}")
message("CXX ID: ${CMAKE_CXX_COMPILER_ID}")


add_subdirectory(src)
