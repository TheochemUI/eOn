# LAMMPS potential: loaded at runtime via dlopen, no compile-time dependency.
# A coverage build links the counter writer with a strong reference.
_lammps_args = _args
if get_option('b_coverage')
    _lammps_args += ['-DEON_GCOV_DUMP']
endif
lammps_pot = library(
    'lammps_pot',
    'LAMMPSPot.cpp',
    'LammpsLoader.cpp',
    dependencies: _deps,
    cpp_args: _lammps_args,
    link_with: _linkto,
    include_directories: _incdirs,
    install: true,
)
