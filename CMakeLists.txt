cmake_minimum_required(VERSION 3.12)
project(TABI-PB VERSION 3.0 LANGUAGES CXX C)


################################################################################
# CMake modules
################################################################################
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})


################################################################################
# Build type
################################################################################
set(default_build_type "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")

  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)

  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release")
endif ()


################################################################################
# Install directories
################################################################################
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)


################################################################################
# OpenACC
################################################################################
option(ENABLE_OPENACC "OpenACC" OFF)

if (ENABLE_OPENACC)
    find_package(OpenACC REQUIRED)
    set(OPENACC_ACCEL_TARGET tesla)
    add_definitions(-DOPENACC_ENABLED)
endif ()


################################################################################
# OpenMP
################################################################################
option(ENABLE_OPENMP "OpenMP" OFF)

if (ENABLE_OPENMP)
    find_package(OpenMP)
    add_definitions(-DOPENMP_ENABLED)
endif ()


################################################################################
# Getting nanoshaper binary
################################################################################
option(GET_NanoShaper "Put NanoShaper executable in the build directory" OFF)

if (GET_NanoShaper)
    set(MESH_ROOT_URL_PATH "https://github.com/lwwilson1/mesh_routines/releases/download/v1.6")
    set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)

    set(NANOSHAPER_FULL_URL_PATH)
    set(NANOSHAPER_RAW_BINARY_NAME "NanoShaper")
    message(STATUS "Copying NanoShaper executable")

    if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")      # Mac OS X specific code
        set(NANOSHAPER_FULL_URL_PATH "${MESH_ROOT_URL_PATH}/NanoShaper_OSX")
    elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")   # Linux specific code
        set(NANOSHAPER_FULL_URL_PATH "${MESH_ROOT_URL_PATH}/NanoShaper_Linux64")
    elseif (${CMAKE_SYSTEM_NAME} MATCHES "Windows") # Windows specific code
        set(NANOSHAPER_FULL_URL_PATH "${MESH_ROOT_URL_PATH}/NanoShaper32.exe")
        set(NANOSHAPER_RAW_BINARY_NAME "NanoShaper.exe")
    endif ()

    if (NOT EXISTS "${EXECUTABLE_OUTPUT_PATH}/${NANOSHAPER_RAW_BINARY_NAME}")
        file(DOWNLOAD ${NANOSHAPER_FULL_URL_PATH} ${EXECUTABLE_OUTPUT_PATH}/${NANOSHAPER_RAW_BINARY_NAME})
        message(STATUS "NanoShaper download: ${NANOSHAPER_FULL_URL_PATH} ${EXECUTABLE_OUTPUT_PATH}/${NANOSHAPER_RAW_BINARY_NAME}")
    endif ()

    install(PROGRAMS ${EXECUTABLE_OUTPUT_PATH}/${NANOSHAPER_RAW_BINARY_NAME} DESTINATION bin)
endif ()


################################################################################
# Setting up src builds
################################################################################
add_subdirectory(src)
