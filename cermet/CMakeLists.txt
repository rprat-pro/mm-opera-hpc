project(OperaHPCBubble)

# behavior law


add_custom_command(
  OUTPUT "src/libBehaviour.so"
  COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/law/Elasticity.mfront ${PROJECT_BINARY_DIR}/law/Norton.mfront ${PROJECT_BINARY_DIR}/law/NortonCr.mfront  ${PROJECT_BINARY_DIR}/law/Mono_UO2_CosH_Jaco.mfront
  DEPENDS ${PROJECT_BINARY_DIR}/law/Elasticity.mfront ${PROJECT_BINARY_DIR}/law/Norton.mfront ${PROJECT_BINARY_DIR}/law/NortonCr.mfront  ${PROJECT_BINARY_DIR}/law/Mono_UO2_CosH_Jaco.mfront
  COMMENT "Performing mfront code generation")

add_custom_target(generate-lib-cermet ALL
  DEPENDS "src/libBehaviour.so")

# build executable

add_executable(cermet-elasticity cermet-elasticity.cxx)
add_executable(cermet-norton cermet-norton.cxx)
add_executable(cermet-norton-bm cermet-norton-bm.cxx)

add_dependencies(cermet-elasticity generate-lib-cermet)
add_dependencies(cermet-norton generate-lib-cermet)
add_dependencies(cermet-norton-bm generate-lib-cermet)

## ELASTICITY
target_link_libraries(cermet-elasticity
 PRIVATE TFELMath TFELException
 PUBLIC mfem-mgis::MFEMMGIS)

target_include_directories(cermet-elasticity
   PUBLIC
   ${TFEL_INCLUDE_PATH} 
   ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_directories(cermet-elasticity
   PUBLIC
   ${TFEL_LIBRARY_PATH}) 

## NORTON
target_link_libraries(cermet-norton
 PRIVATE TFELMath TFELException
 PUBLIC mfem-mgis::MFEMMGIS)

target_include_directories(cermet-norton
   PUBLIC
   ${TFEL_INCLUDE_PATH} 
   ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_directories(cermet-norton
   PUBLIC
   ${TFEL_LIBRARY_PATH}) 

## NORTON BRUNO MICHEL
target_link_libraries(cermet-norton-bm
 PRIVATE TFELMath TFELException
 PUBLIC mfem-mgis::MFEMMGIS)

target_include_directories(cermet-norton-bm
   PUBLIC
   ${TFEL_INCLUDE_PATH} 
   ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_directories(cermet-norton-bm
   PUBLIC
   ${TFEL_LIBRARY_PATH}) 

# copy and install

file(COPY law mesh DESTINATION ${PROJECT_BINARY_DIR})

add_test(NAME CermetElasticity COMMAND cermet-elasticity --mesh ${PROJECT_BINARY_DIR}/mesh/cermet-mini.msh)
add_test(NAME CermetNorton COMMAND cermet-norton --mesh ${PROJECT_BINARY_DIR}/mesh/cermet-mini.msh)

install(TARGETS cermet-elasticity DESTINATION bin)
install(TARGETS cermet-norton DESTINATION bin)
