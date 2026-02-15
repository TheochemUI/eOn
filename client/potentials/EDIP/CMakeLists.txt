add_library(
  EDIP SHARED
  EDIP.cpp
  edipFortran.f90)
target_link_libraries(EDIP PRIVATE eoncbase ${EON_DEPS})
target_include_directories(EDIP PRIVATE ${EON_INCLUDE_DIRS})
if(EON_FORTRAN_FLAGS)
  target_compile_options(EDIP PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${EON_FORTRAN_FLAGS}>)
endif()
install(TARGETS EDIP)
