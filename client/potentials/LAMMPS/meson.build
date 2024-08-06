libdir = meson.current_source_dir()
lammps_so = cppc.find_library('lammps', dirs: libdir)  # ./liblammps.so
potss.add(when: ['LAMMPS_POT', lammps_so], if_true: files('LAMMPSPot.cpp'))
