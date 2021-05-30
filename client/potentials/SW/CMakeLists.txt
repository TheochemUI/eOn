add_library(
  SW SHARED
  SW.cpp
  SWFortran.f90
  )

# Install
install(TARGETS SW)
