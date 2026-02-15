add_library(
  ASE_NWCHEM SHARED
  ASE_NWCHEM.cpp
  ../../EnvHelpers.cc
  ../../fpe_handler.cpp)
target_link_libraries(ASE_NWCHEM PRIVATE eoncbase ${EON_DEPS})
target_include_directories(ASE_NWCHEM PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS ASE_NWCHEM)
