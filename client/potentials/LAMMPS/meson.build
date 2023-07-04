libdir = meson.current_source_dir()
lammps_so = cppc.find_library('lammps', dirs : libdir) # ./liblammps.so
_deps += [ lammps_so ]
lammps_pot = shared_library('lammps_pot',
           'LAMMPS.cpp',
           dependencies : _deps,
           cpp_args : _args,
           install : true)
