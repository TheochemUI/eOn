# LAMMPS potential: loaded at runtime via dlopen, no compile-time dependency.
# -ldl is added globally to _deps in client/meson.build (file-scope dl_dep)
# so every loader-using target links it without per-potential boilerplate.
lammps_pot = library('lammps_pot',
    'LAMMPSPot.cpp',
    'LammpsLoader.cpp',
    'LammpsBundle.cpp',
    dependencies: _deps,
    cpp_args: _args,
    link_with: _linkto,
    include_directories: _incdirs,
    install: true)
