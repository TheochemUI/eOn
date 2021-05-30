add_library(
  Lenosky SHARED
  Lenosky.cpp
  lenoskyFortran.f90
  )

# Install
install(TARGETS Lenosky)
