# LAMMPS potential: liblammps is loaded at runtime via dlopen.
# A coverage build links the counter writer with a strong reference.
# This plugin itself is linked into eonclib. On Windows a shared library
# with no dllexport produces no import library.
_lammps_args = _args
if get_option('b_coverage')
    _lammps_args += ['-DEON_GCOV_DUMP']
endif
_lammps_opts = []
_lammps_install = true
if is_windows
    _lammps_opts = ['default_library=static']
    _lammps_install = false
endif
lammps_pot = library(
    'lammps_pot',
    'LAMMPSPot.cpp',
    'LammpsLoader.cpp',
    dependencies: _deps,
    cpp_args: _lammps_args,
    link_with: _linkto,
    include_directories: _incdirs,
    override_options: _lammps_opts,
    install: _lammps_install,
)
