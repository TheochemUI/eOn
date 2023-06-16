libdir = meson.current_source_dir()
lammps_so = cppc.find_library('lammps', dirs : libdir) # ./liblammps.so
lammps_pot = shared_library('lammps_pot',
           'LAMMPS.cpp',
           dependencies : [eon_deps, lammps_so],
           cpp_args : eon_extra_args,
           install : true)
