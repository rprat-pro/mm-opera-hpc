project(OperaHPCBubble)

# behavior law


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
 PRIVATE TFELMath TFELException
 PUBLIC mfem-mgis::MFEMMGIS)

target_include_directories(test-bubble
   PUBLIC
   ${TFEL_INCLUDE_PATH} 
   ${CMAKE_CURRENT_SOURCE_DIR}/include
   ${HDF5_INCLUDE_DIRS})

target_link_directories(test-bubble
   PUBLIC
   ${TFEL_LIBRARY_PATH}
   ${HDF5_LIBRARIES}) 

if(TARGET HDF5::hdf5_cpp AND TARGET HDF5::hdf5)
    message(STATUS "Linking with HDF5 imported targets: HDF5::hdf5_cpp and HDF5::hdf5")
    # Order usually doesn't strictly matter with imported targets if dependencies are set correctly,
    # but listing the C++ one first or only the C++ one (if it pulls in the C one) is common.
    target_link_libraries(test-bubble PRIVATE HDF5::hdf5_cpp HDF5::hdf5)
    # Or, sometimes HDF5::hdf5_cpp itself will link HDF5::hdf5, so just:
    # target_link_libraries(test-bubble PRIVATE HDF5::hdf5_cpp)
    # Or sometimes a single HDF5::HDF5 target is provided that covers all found components:
    # target_link_libraries(test-bubble PRIVATE HDF5::HDF5)
elseif(HDF5_FOUND AND DEFINED HDF5_CXX_LIBRARIES AND DEFINED HDF5_C_LIBRARIES)
    # Fallback to using the older variable-based approach if imported targets aren't available
    # or if using an older FindHDF5.cmake module from CMake itself.
    message(STATUS "Linking with HDF5 library variables:")
    message(STATUS "  HDF5 CXX Libraries: ${HDF5_CXX_LIBRARIES}")
    message(STATUS "  HDF5 C Libraries: ${HDF5_C_LIBRARIES}")
    message(STATUS "  HDF5 Include Dirs: ${HDF5_INCLUDE_DIRS}")

    # Include HDF5 header directories
    target_include_directories(test-bubble PRIVATE ${HDF5_INCLUDE_DIRS})

    # Link the libraries. The order can be important here: C++ library first, then C library.
    # ${HDF5_CXX_LIBRARIES} should contain libhdf5_cpp.
    # ${HDF5_C_LIBRARIES} should contain libhdf5.
    target_link_libraries(test-bubble PRIVATE ${HDF5_CXX_LIBRARIES} ${HDF5_C_LIBRARIES})
else()
    message(FATAL_ERROR "HDF5 C++ libraries not found or not configured correctly. "
                        "Please ensure HDF5 (with C++ support) is installed and findable by CMake. "
                        "HDF5_FOUND=${HDF5_FOUND}, HDF5_CXX_LIBRARIES=${HDF5_CXX_LIBRARIES}")
endif()


# copy and install

file(COPY Elasticity.mfront mesh DESTINATION ${PROJECT_BINARY_DIR})

add_test(NAME BubbleTest COMMAND test-bubble --mesh mesh/mesh_sphere.msh --bubble-file mesh/single_bubble_ci.txt)

install(TARGETS test-bubble DESTINATION bin)
