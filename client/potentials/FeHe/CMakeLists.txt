add_library(
  FeHe SHARED
  FeHe.cpp
  feforce.f
  )

# Install
install(TARGETS FeHe)
