set(basePotentials
    Aluminum
    EAM
    EDIP
    EMT
    ExtPot
    FeHe
    GPRPotential
    # LAMMPS # Needs external library
    Lenosky
    LJ
    LJCluster
    Morse
    NewPot
    VASP
    # PyAMFF # Needs external library
    QSC
    SW
    Tersoff
    AMS
    AMS_IO
    Water
    Water_H
    Water_Pt)

if(${WIN32})
  # No VASP for Windows
  string(REPLACE "VASP" "" basePotentials "${basePotentials}")
endif()

if(${MPI})
  # Add when MPI is avaliable, ; is important
  string(APPEND basePotentials ";MPIPot")
endif()

if(${NO_FORTRAN})
  # Kill all Fortran based potentials
  string(REPLACE "Aluminum" "" basePotentials "${basePotentials}")
  string(REPLACE "Lenosky" "" basePotentials "${basePotentials}")
  string(REPLACE "SW" "" basePotentials "${basePotentials}")
  string(REPLACE "Tersoff" "" basePotentials "${basePotentials}")
  string(REPLACE "EDIP" "" basePotentials "${basePotentials}")
  string(REPLACE "Water_H" "" basePotentials "${basePotentials}")
  string(REPLACE "FeHe" "" basePotentials "${basePotentials}")
endif()

foreach(pot ${basePotentials})
  add_subdirectory(${pot})
endforeach()

add_library(potentials INTERFACE)
target_link_libraries(potentials INTERFACE ${basePotentials})
