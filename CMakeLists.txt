cmake_minimum_required(VERSION 3.22)
project(opera-hpc
        VERSION 0.1.0
        DESCRIPTION "Evolution of a population of bubble in a brittle material")

# ===========================
# === CMake customization ===
# ===========================
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Coverage RelWithDebInfo Debug)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)
include(cmake/modules/tfel.cmake)

find_package(MFEMMGIS REQUIRED)

enable_testing()

add_subdirectory(polycrystal)
add_subdirectory(bubble)
add_subdirectory(cermet)

