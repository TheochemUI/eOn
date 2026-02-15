add_library(
  Lenosky SHARED
  Lenosky.cpp
  lenoskyFortran.f90)
target_link_libraries(Lenosky PRIVATE eoncbase ${EON_DEPS})
target_include_directories(Lenosky PRIVATE ${EON_INCLUDE_DIRS})
if(EON_FORTRAN_FLAGS)
  target_compile_options(Lenosky PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${EON_FORTRAN_FLAGS}>)
endif()
install(TARGETS Lenosky)
