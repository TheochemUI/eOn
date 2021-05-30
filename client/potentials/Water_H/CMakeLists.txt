add_library(
  Water_H SHARED
  potH_H2O.f90
  Tip4p_H.cpp
  )

# Install
install(TARGETS Water_H)
