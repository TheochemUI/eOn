# LAMMPS potential: liblammps is loaded at runtime via dlopen.
# A coverage build links the counter writer with a strong reference.
# This plugin is linked into eonclib. library() follows default_library,
# so a Windows shared build emits a DLL with no import library.
_lammps_args = _args
if get_option('b_coverage')
    _lammps_args += ['-DEON_GCOV_DUMP']
endif
if is_windows
    lammps_pot = static_library(
        'lammps_pot',
        'LAMMPSPot.cpp',
        'LammpsLoader.cpp',
        dependencies: _deps,
        cpp_args: _lammps_args,
        link_with: _linkto,
        include_directories: _incdirs,
        install: get_option('default_library') != 'shared',
    )
else
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
endif
