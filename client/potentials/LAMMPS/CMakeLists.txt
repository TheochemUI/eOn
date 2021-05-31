add_library(
  LAMMPS SHARED
  LAMMPS.cpp
  )
# Real LAMMPS (must be provided)
add_library(liblammps SHARED IMPORTED)
set_property(TARGET liblammps PROPERTY IMPORTED_LOCATION "./liblammps.so")
target_link_libraries(LAMMPS liblammps)

# Install
install(TARGETS LAMMPS)
