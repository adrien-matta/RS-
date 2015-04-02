add_library(RSCore SHARED RSDifferentialCrossSection.cxx RSShellModelState.cxx RSShellModelCollection.cxx RSExperimentalCrossSection.cxx RSOverlapFunction.cxx)
target_link_libraries(RSCore ${ROOT_LIBRARIES})
