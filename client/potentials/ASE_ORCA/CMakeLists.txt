add_library(
  ASE_ORCA SHARED
  ASE_ORCA.cpp
  ../../EnvHelpers.cc
  ../../fpe_handler.cpp)
target_link_libraries(ASE_ORCA PRIVATE eoncbase ${EON_DEPS})
target_include_directories(ASE_ORCA PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS ASE_ORCA)
