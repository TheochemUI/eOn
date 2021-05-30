add_library(
  LAMMPS SHARED
  LAMMPS.cpp
  )

# Install
install(TARGETS LAMMPS)
