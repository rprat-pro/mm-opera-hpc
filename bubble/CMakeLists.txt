project(OperaHPCBubble)

# behavior law
file(COPY verification/mmm_vs_analytical.py DESTINATION ${PROJECT_BINARY_DIR})


add_custom_command(
  OUTPUT "src/libBehaviour.so"
  COMMAND mfront --obuild --interface=generic ${PROJECT_BINARY_DIR}/Elasticity.mfront
  DEPENDS ${PROJECT_BINARY_DIR}/Elasticity.mfront
  COMMENT "Performing mfront code generation")

add_custom_target(lib-elasticity ALL
  DEPENDS "src/libBehaviour.so")

# build executable

file( GLOB_RECURSE SRCS *.cxx )

add_executable(test-bubble
  ${SRCS}
)

add_dependencies(test-bubble lib-elasticity)

target_link_libraries(test-bubble
 PRIVATE tfel::TFELMath tfel::TFELException
 PUBLIC mfem-mgis::MFEMMGIS)

target_include_directories(test-bubble
   PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# copy and install

file(COPY Elasticity.mfront mesh DESTINATION ${PROJECT_BINARY_DIR})

add_test(NAME BubbleTest COMMAND test-bubble --mesh mesh/mesh_sphere.msh --bubble-file mesh/single_bubble_ci.txt)


install(TARGETS test-bubble DESTINATION bin)
