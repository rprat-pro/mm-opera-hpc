add_library(MMOPERAHPC
        GrainOrientations.cxx
        MacroscropicElasticMaterialProperties.cxx
        UniaxialMacroscopicStressPeriodicSimulation.cxx)

target_include_directories(MMOPERAHPC
   PUBLIC
   ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(MMOPERAHPC
 PUBLIC mfem-mgis::MFEMMGIS)
