set(basePotentials
    EAM
    EMT
    ExtPot
    LJ
    LJCluster
    Morse
    NewPot
    VASP
    )

if(${WIN32})
  # No VASP for Windows
  string(REPLACE "VASP" "" basePotentials "${basePotentials}")
endif()

if(${MPI})
  # Add when MPI is avaliable, ; is important
  string(APPEND basePotentials ";MPIPot")
endif()

if(${WITH_AMS})
  # Add when AMS is avaliable, ; is important
  string(APPEND basePotentials ";AMS;AMS_IO")
endif()

if(${WITH_WATER})
  string(APPEND basePotentials ";Water;Water_H;Water_Pt")
endif()

if(${WITH_LAMMPS})
  string(APPEND basePotentials ";LAMMPS")
endif()

if(${WITH_FORTRAN})
  string(APPEND basePotentials ";Aluminum;Water_H;Lenosky;SW;Tersoff;EDIP;FeHe;SW")
endif()

if(${WITH_PYTHON})
  string(APPEND basePotentials ";PyAMFF;QSC")
endif()

if(${WITH_GPRD})
  string(APPEND basePotentials ";GPRPotential")
endif()

foreach(pot ${basePotentials})
  add_subdirectory(${pot})
endforeach()

add_library(potentials INTERFACE)
target_link_libraries(potentials INTERFACE ${basePotentials})

# Install
install(TARGETS potentials)
