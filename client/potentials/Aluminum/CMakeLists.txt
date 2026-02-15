add_library(
  Aluminum SHARED
  dfrhoDblexp.f
  embedenergy.f
  fofrhoDblexp.f
  forces.f
  gagafeDblexp.f
  potinit.f
  sumembforce.f
  Aluminum.cpp)
target_link_libraries(Aluminum PRIVATE eoncbase ${EON_DEPS})
target_include_directories(Aluminum PRIVATE ${EON_INCLUDE_DIRS})
if(EON_FORTRAN_FLAGS)
  target_compile_options(Aluminum PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${EON_FORTRAN_FLAGS}>)
endif()
install(TARGETS Aluminum)
