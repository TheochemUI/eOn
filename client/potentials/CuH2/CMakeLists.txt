add_library(
  CuH2 SHARED
  CuH2.cpp
  eam_isoc.f90
  eamroutines.f90)
target_link_libraries(CuH2 PRIVATE eoncbase ${EON_DEPS})
target_include_directories(CuH2 PRIVATE ${EON_INCLUDE_DIRS})
if(EON_FORTRAN_FLAGS)
  target_compile_options(CuH2 PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${EON_FORTRAN_FLAGS}>)
endif()
install(TARGETS CuH2)
