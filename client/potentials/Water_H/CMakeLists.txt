add_library(
  Water_H SHARED
  potH_H2O.f90
  Tip4p_H.cpp)
target_link_libraries(Water_H PRIVATE eoncbase ${EON_DEPS})
target_include_directories(Water_H PRIVATE ${EON_INCLUDE_DIRS})
if(EON_FORTRAN_FLAGS)
  target_compile_options(Water_H PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${EON_FORTRAN_FLAGS}>)
endif()
install(TARGETS Water_H)
