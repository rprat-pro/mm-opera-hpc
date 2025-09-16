project(OperaHPCBubble)

# behavior law
add_custom_command(
  OUTPUT "src/libBehaviour.so"
  COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/law/NortonCr.mfront  ${PROJECT_BINARY_DIR}/law/Mono_UO2_CosH_Jaco.mfront
  DEPENDS ${PROJECT_BINARY_DIR}/law/NortonCr.mfront  ${PROJECT_BINARY_DIR}/law/Mono_UO2_CosH_Jaco.mfront
  COMMENT "Performing mfront code generation")

add_custom_target(generate-lib-cermet ALL
  DEPENDS "src/libBehaviour.so")

# build executable
add_executable(cermet cermet.cxx)
add_dependencies(cermet generate-lib-cermet)

target_link_libraries(cermet
 PUBLIC MMOPERAHPC)

# copy and install

file(COPY law mesh DESTINATION ${PROJECT_BINARY_DIR})
file(COPY mesh/vectors_5grains.txt DESTINATION ${PROJECT_BINARY_DIR})

add_test(NAME CermetNullStrainXXYY COMMAND cermet --duration 0.5 --nstep 1 --bcs_type 0)
add_test(NAME CermetImposeGradXXYY COMMAND cermet --duration 0.5 --nstep 1 --bcs_type 1)

install(TARGETS cermet DESTINATION bin)
