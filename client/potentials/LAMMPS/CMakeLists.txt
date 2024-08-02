# Real LAMMPS (must be provided)
add_library(truelmp SHARED IMPORTED)
set_property(TARGET truelmp PROPERTY IMPORTED_LOCATION "./liblammps.so")

# Wrapper
add_library(
  LAMMPS SHARED
  LAMMPS.cpp
  )

# Link
target_link_libraries(LAMMPS truelmp)

# Install
install(TARGETS LAMMPS)
