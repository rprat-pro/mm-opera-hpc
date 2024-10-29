cmake_minimum_required(VERSION 3.22)
project(opera-hpc
        VERSION 0.1.0
        DESCRIPTION "Evolution of a population of bubble in a brittle material")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include(cmake/modules/tfel.cmake)

find_package(MFEM REQUIRED)
find_package(MFrontGenericInterface REQUIRED)
find_package(MFEMMGIS REQUIRED)
#find_package(TFELMath REQUIRED)
#include_directories(/home/rp269144/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-11.4.0/tfel-master-u6q5krrbckkhn3ahzrwjils5oijqvfzv/include/)

enable_testing()

add_subdirectory(bubble)
add_subdirectory(polycrystal)

