add_library(
  EMT SHARED
  EffectiveMediumTheory.cpp
  Asap/Atoms.cpp
  Asap/NeighborList.cpp
  Asap/CellList.cpp
  Asap/EMT.cpp
  Asap/SuperCell.cpp
  Asap/EMTDefaultParameterProvider.cpp
  Asap/Vec.cpp
  Asap/Exception.cpp
  Asap/GhostAtoms.cpp
  Asap/EMTRasmussenParameterProvider.cpp)

# Install
install(TARGETS EMT)
