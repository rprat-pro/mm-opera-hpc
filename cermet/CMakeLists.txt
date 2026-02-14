project(OperaHPCCermetExample)

# Copy
file(COPY law mesh results DESTINATION ${PROJECT_BINARY_DIR})
file(COPY mesh/vectors_5grains.txt check_cermet_restults.py plot_cermet_results.py DESTINATION ${PROJECT_BINARY_DIR})

# behavior law
add_custom_command(
  OUTPUT "src/libBehaviour.so"
  COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/law/NortonCr.mfront
          ${PROJECT_BINARY_DIR}/law/Mono_UO2_CosH_Jaco.mfront
          ${PROJECT_BINARY_DIR}/law/Mono_UO2_CosH_Jaco3.mfront
  DEPENDS ${PROJECT_BINARY_DIR}/law/NortonCr.mfront
          ${PROJECT_BINARY_DIR}/law/Mono_UO2_CosH_Jaco.mfront
          ${PROJECT_BINARY_DIR}/law/Mono_UO2_CosH_Jaco3.mfront
  COMMENT "Performing mfront code generation")

add_custom_target(generate-lib-cermet ALL
  DEPENDS "src/libBehaviour.so")

## Build executable
add_executable(cermet cermet.cxx)
target_link_libraries(cermet
 PUBLIC MMOPERAHPC)
add_dependencies(cermet generate-lib-cermet)

## Test
add_dependencies(check cermet generate-lib-cermet)
add_test(NAME CermetNullStrainXXYY COMMAND cermet --duration 0.4 --nstep 1 --postprocessing 0)
add_test(NAME CermetImposeGradXXYY COMMAND cermet --duration 0.4 --nstep 1 --postprocessing 0)

## Install
install(TARGETS cermet DESTINATION bin)
