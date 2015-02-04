add_library(RSReaction SHARED RSDifferentialCrossSection.h RSDifferentialCrossSection.cxx )
target_link_libraries(RSReaction ${ROOT_LIBRARIES})
