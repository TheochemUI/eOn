# LAMMPS potential: loaded at runtime via dlopen, no compile-time dependency
lammps_pot = library('lammps_pot',
    'LAMMPSPot.cpp',
    'LammpsLoader.cpp',
    dependencies: _deps,
    cpp_args: _args,
    link_with: _linkto,
    include_directories: _incdirs,
    install: true)
