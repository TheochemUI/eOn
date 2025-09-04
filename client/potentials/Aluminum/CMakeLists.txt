add_library(
  Aluminum SHARED
  dfrhoDblexp.f
  embedenergy.f
  fofrhoDblexp.f
  forces.f
  gagafeDblexp.f
  potinit.f
  sumembforce.f
  Aluminum.cpp
  )

# Install
install(TARGETS Aluminum)
