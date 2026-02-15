add_library(
  SW SHARED
  SW.cpp
  SWFortran.f90)
target_link_libraries(SW PRIVATE eoncbase ${EON_DEPS})
target_include_directories(SW PRIVATE ${EON_INCLUDE_DIRS})
if(EON_FORTRAN_FLAGS)
  target_compile_options(SW PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${EON_FORTRAN_FLAGS}>)
endif()
install(TARGETS SW)
