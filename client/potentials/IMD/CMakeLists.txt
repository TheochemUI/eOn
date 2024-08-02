add_library(
  IMD SHARED
  IMD.cpp
  IMD.h
  filesForIMD/imd_cg.c
  )

# Install
install(TARGETS IMD)
