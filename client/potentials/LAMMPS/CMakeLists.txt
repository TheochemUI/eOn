find_library(LAMMPS_LIB lammps QUIET)

add_library(LAMMPS SHARED LAMMPSPot.cpp)
target_include_directories(LAMMPS PRIVATE ${EON_INCLUDE_DIRS})

if(LAMMPS_LIB)
  target_link_libraries(LAMMPS PRIVATE eoncbase ${EON_DEPS} ${LAMMPS_LIB})
else()
  # Fallback: look for liblammps.so in the current source directory
  add_library(truelmp SHARED IMPORTED)
  set_property(TARGET truelmp PROPERTY IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/liblammps.so")
  target_link_libraries(LAMMPS PRIVATE eoncbase ${EON_DEPS} truelmp)
endif()

install(TARGETS LAMMPS)
