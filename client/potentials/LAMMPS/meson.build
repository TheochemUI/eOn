libdir = meson.current_source_dir()
lammps_so = cppc.find_library('lammps', dirs : libdir) # ./liblammps.so
lammps_pot = shared_library('lammps_pot',
           'LAMMPSPot.cpp',
           dependencies : [ _deps, lammps_so],
           cpp_args : _args,
           link_with : _linkto,
           include_directories: _incdirs,
           install : true)
