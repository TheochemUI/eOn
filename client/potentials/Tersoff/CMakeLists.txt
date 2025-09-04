add_library(
  Tersoff SHARED
  Tersoff.cpp
  tersoffFortran.f90
  )

# Install
install(TARGETS Tersoff)
