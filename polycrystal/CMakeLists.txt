project(UniaxialPolycrystalExample)

# Copy needed file for execution
file(COPY
  ${PROJECT_SOURCE_DIR}/law/Mono_UO2_CosH_Jaco.mfront
  ${PROJECT_SOURCE_DIR}/law/Mono_UO2_CosH_Jaco3.mfront
  DESTINATION ${PROJECT_BINARY_DIR})

file(COPY law mesh results plot_polycrystal_results.py check_polycrystal_restults.py DESTINATION ${PROJECT_BINARY_DIR})
  
  # Use mfront to generate libBehaviour.so
add_custom_command(
  OUTPUT "src/libBehaviour.so"
  COMMAND mfront --obuild --interface=generic
          ${PROJECT_BINARY_DIR}/law/Mono_UO2_CosH_Jaco.mfront
          ${PROJECT_BINARY_DIR}/law/Mono_UO2_CosH_Jaco3.mfront
  DEPENDS ${PROJECT_BINARY_DIR}/law/Mono_UO2_CosH_Jaco.mfront
          ${PROJECT_BINARY_DIR}/law/Mono_UO2_CosH_Jaco3.mfront
  COMMENT "Performing mfront code generation")
  
add_custom_target(generate-libBehaviour-polycrystal ALL
  DEPENDS "src/libBehaviour.so")

add_executable(uniaxial-polycrystal UniaxialPolycrystal.cxx)
target_link_libraries(uniaxial-polycrystal MMOPERAHPC)
target_include_directories(uniaxial-polycrystal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
    
# testing
add_dependencies(check uniaxial-polycrystal generate-libBehaviour-polycrystal)
add_test(NAME UniaxialPolycrystalTest COMMAND uniaxial-polycrystal --disable-post-processing --nstep 1 --duration 0.1)
add_test(NAME UniaxialPolycrystalTestPar COMMAND mpirun -n 2 uniaxial-polycrystal --disable-post-processing --nstep 1 --duration 0.1)
add_test(NAME PolycrystalTestWithPostProPar COMMAND mpirun -n 2 uniaxial-polycrystal --enable-export-von-Mises-stress --enable-export-first_eigen_stress --nstep 1 --duration 0.01)
