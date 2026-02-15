lammps_so = dependency('liblammps', required: false)
if not lammps_so.found()
    libdir = meson.current_source_dir()
    lammps_so = cppc.find_library('lammps', dirs: libdir)
endif
lammps_pot = library('lammps_pot',
           'LAMMPSPot.cpp',
           dependencies : [ _deps, lammps_so],
           cpp_args : _args,
           link_with : _linkto,
           include_directories: _incdirs,
           install : true)
