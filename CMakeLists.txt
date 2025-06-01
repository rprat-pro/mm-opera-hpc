cmake_minimum_required(VERSION 3.22)
project(opera-hpc
        VERSION 0.1.0
        DESCRIPTION "Evolution of a population of bubble in a brittle material")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

include(cmake/modules/tfel.cmake)

find_package(MFEMMGIS REQUIRED)
#find_package(TFELMath REQUIRED)
#include_directories(/home/rp269144/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-11.4.0/tfel-master-u6q5krrbckkhn3ahzrwjils5oijqvfzv/include/)

enable_testing()

add_subdirectory(bubble)
add_subdirectory(polycrystal)

