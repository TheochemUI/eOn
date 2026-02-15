add_library(
  FeHe SHARED
  FeHe.cpp
  feforce.f)
target_link_libraries(FeHe PRIVATE eoncbase ${EON_DEPS})
target_include_directories(FeHe PRIVATE ${EON_INCLUDE_DIRS})
if(EON_FORTRAN_FLAGS)
  target_compile_options(FeHe PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${EON_FORTRAN_FLAGS}>)
endif()
install(TARGETS FeHe)
