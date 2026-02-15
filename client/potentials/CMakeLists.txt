# Base potentials (always built)
set(basePotentials
  EAM
  EMT
  ExtPot
  IMD
  LJ
  LJCluster
  Morse
  ZBL
  NewPot)

if(NOT WIN32)
  list(APPEND basePotentials SocketNWChem)
endif()

# ---------------------- Conditional potentials

if(WITH_VASP AND NOT WIN32)
  list(APPEND basePotentials VASP)
endif()

if(WITH_MPI)
  list(APPEND basePotentials MPIPot)
endif()

if(WITH_AMS)
  list(APPEND basePotentials AMS AMS_IO)
endif()

if(WITH_WATER)
  list(APPEND basePotentials Water Water_Pt)
endif()

if(WITH_LAMMPS)
  list(APPEND basePotentials LAMMPS)
endif()

if(WITH_FORTRAN)
  list(APPEND basePotentials Aluminum EDIP FeHe Lenosky SW Tersoff)
  if(WITH_WATER)
    list(APPEND basePotentials Water_H)
  endif()
endif()

if(WITH_CUH2)
  list(APPEND basePotentials CuH2)
endif()

if(WITH_PYTHON OR EON_NEED_PYTHON_EMBED)
  list(APPEND basePotentials PyAMFF QSC)
endif()

if(WITH_GPRD)
  list(APPEND basePotentials GPRPotential)
endif()

if(WITH_XTB)
  list(APPEND basePotentials XTBPot)
endif()

if(WITH_ASE AND EON_NEED_PYTHON_EMBED)
  list(APPEND basePotentials ASE)
endif()

if(WITH_ASE_ORCA AND EON_NEED_PYTHON_EMBED)
  list(APPEND basePotentials ASE_ORCA)
endif()

if(WITH_ASE_NWCHEM AND EON_NEED_PYTHON_EMBED)
  list(APPEND basePotentials ASE_NWCHEM)
endif()

if(WITH_CATLEARN AND EON_NEED_PYTHON_EMBED)
  list(APPEND basePotentials CatLearnPot)
endif()

if(WITH_METATOMIC)
  list(APPEND basePotentials Metatomic)
endif()

# ---------------------- Build all selected potentials
foreach(pot ${basePotentials})
  add_subdirectory(${pot})
endforeach()

# Interface library linking all potentials
add_library(potentials INTERFACE)
target_link_libraries(potentials INTERFACE ${basePotentials})

install(TARGETS potentials)
