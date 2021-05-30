add_library(
  EDIP SHARED
  EDIP.cpp
  edipFortran.f90
  )

# Install
install(TARGETS EDIP)
