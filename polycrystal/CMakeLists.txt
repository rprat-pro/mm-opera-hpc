cmake_minimum_required(VERSION 3.13)
project(UniaxialPolycristalExample)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CTest)
# testing
set(CTEST_CONFIGURATION_TYPE "${JOB_BUILD_CONFIGURATION}")
# (must be placed *before* any add_subdirectory, cmake bug ?)
enable_testing()
if(CMAKE_CONFIGURATION_TYPES)
   if (NOT TARGET check)
      add_custom_target(check COMMAND 
         ${CMAKE_CTEST_COMMAND} -T test -C $<CONFIGURATION>)
   endif(NOT TARGET check)
else(CMAKE_CONFIGURATION_TYPES)
   if (NOT TARGET check)
      add_custom_target(check COMMAND 
      ${CMAKE_CTEST_COMMAND} -T test )
   endif(NOT TARGET check)
endif(CMAKE_CONFIGURATION_TYPES)


# Define a test that will be called during the "make check"

# Copy needed file for execution
file(COPY
  ${PROJECT_SOURCE_DIR}/Mono_UO2_CosH_Jaco.mfront
  ${PROJECT_SOURCE_DIR}/Elasticity.mfront
  ${PROJECT_SOURCE_DIR}/ElasticityOrtho.mfront
  ${PROJECT_SOURCE_DIR}/SaintVenantKirchhoffElasticity.mfront
  ${PROJECT_SOURCE_DIR}/periodic-cube.msh
  ${PROJECT_SOURCE_DIR}/cube_2mat_per.mesh
  ${PROJECT_SOURCE_DIR}/n8-id1.msh
  ${PROJECT_SOURCE_DIR}/n8-id1-vecteurs.txt
  DESTINATION ${PROJECT_BINARY_DIR})
  
  # Use mfront to generate libBehaviour.so
add_custom_command(
  OUTPUT "src/libBehaviour.so"
  COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/Mono_UO2_CosH_Jaco.mfront ${PROJECT_BINARY_DIR}/Elasticity.mfront ${PROJECT_BINARY_DIR}/SaintVenantKirchhoffElasticity.mfront ${PROJECT_BINARY_DIR}/ElasticityOrtho.mfront
  DEPENDS ${PROJECT_BINARY_DIR}/Mono_UO2_CosH_Jaco.mfront ${PROJECT_BINARY_DIR}/Elasticity.mfront ${PROJECT_BINARY_DIR}/SaintVenantKirchhoffElasticity.mfront ${PROJECT_BINARY_DIR}/ElasticityOrtho.mfront
  COMMENT "Performing mfront code generation")
  
add_custom_target(generate-libBehaviour-polycrystal ALL
  DEPENDS "src/libBehaviour.so")

function(build_polycrystal binFile srcFile)
	add_executable(${binFile} ${srcFile})
	target_link_libraries(${binFile} mfem-mgis::MFEMMGIS)
	target_include_directories(${binFile}
  	SYSTEM
  	PUBLIC "$<BUILD_INTERFACE:${MFEM_INCLUDE_DIRS}>"
  	PUBLIC "$<INSTALL_INTERFACE:${MFEM_INCLUDE_DIRS}>")
	target_include_directories(${binFile} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
endfunction(build_polycrystal)

build_polycrystal(uniaxial-polycristal UniaxialPolycristal.cxx)
build_polycrystal(uniaxial-elastic UniaxialElastic.cxx)
build_polycrystal(uniaxial-nonlinearelastic UniaxialNonLinearElastic.cxx)
build_polycrystal(uniaxial-elastic-ortho UniaxialElasticOrtho.cxx)
    
add_dependencies(check uniaxial-polycristal uniaxial-elastic uniaxial-elastic-ortho uniaxial-nonlinearelastic generate-libBehaviour-polycrystal)

set(test "UniaxialPolycristalTest")
add_test(NAME ${test} COMMAND uniaxial-polycristal --postprocessings 0)

set(test "UniaxialElasticTest")
add_test(NAME ${test} COMMAND uniaxial-elastic --postprocessings 0)

set(test "UniaxialNonLinearElasticTest")
add_test(NAME ${test} COMMAND uniaxial-nonlinearelastic --postprocessings 0)

