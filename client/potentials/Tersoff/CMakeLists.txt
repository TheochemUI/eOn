add_library(
  Tersoff SHARED
  Tersoff.cpp
  tersoffFortran.f90)
target_link_libraries(Tersoff PRIVATE eoncbase ${EON_DEPS})
target_include_directories(Tersoff PRIVATE ${EON_INCLUDE_DIRS})
if(EON_FORTRAN_FLAGS)
  target_compile_options(Tersoff PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${EON_FORTRAN_FLAGS}>)
endif()
install(TARGETS Tersoff)
