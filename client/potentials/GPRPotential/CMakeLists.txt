add_library(GPRPotential SHARED GPRPotential.cpp)
target_link_libraries(GPRPotential PRIVATE eoncbase ${EON_DEPS})
target_include_directories(GPRPotential PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS GPRPotential)
