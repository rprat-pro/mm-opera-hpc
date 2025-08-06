project(OperaHPCBubble)

# behavior law


add_custom_command(
  OUTPUT "src/libBehaviour.so"
  COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/law/Elasticity.mfront
  DEPENDS ${PROJECT_BINARY_DIR}/law/Elasticity.mfront
  COMMENT "Performing mfront code generation")

add_custom_target(lib-elasticity ALL
  DEPENDS "src/libBehaviour.so")

# build executable

file( GLOB_RECURSE SRCS *.cxx )

add_executable(cermet-elasticity
  ${SRCS}
)

add_dependencies(cermet-elasticity lib-elasticity)

target_link_libraries(cermet-elasticity
 PRIVATE TFELMath TFELException
 PUBLIC mfem-mgis::MFEMMGIS)

target_include_directories(test-bubble
   PUBLIC
   ${TFEL_INCLUDE_PATH} 
   ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_directories(test-bubble
   PUBLIC
   ${TFEL_LIBRARY_PATH}) 


# copy and install

file(COPY Elasticity.mfront mesh DESTINATION ${PROJECT_BINARY_DIR})

add_test(NAME CermetElasticity COMMAND cermet-elasticity --mesh ${PROJECT_BINARY_DIR}/mesh/cermet-mini.msh)

install(TARGETS cermet-elasticity DESTINATION bin)
