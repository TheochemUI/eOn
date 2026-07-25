# Fortran-only shared module: loaded at runtime via dlopen/LoadLibrary.
# Neighbour list built via the vendored vesin Fortran interface.
eon_sw = shared_module(
    'eon_sw',
           'SWFortran.f90',
    fortran_args: _fargs,
    dependencies: [vesin_f_dep],
    link_args: _fortran_pot_link_args
    + (
        is_windows ? [
            '-Wl,/DEF:'
            + join_paths(
                meson.project_source_root(),
                'client',
                'potentials',
                'win_exports',
                'eon_sw.def',
            ),
        ] : []
    ),
    install: true,
)

# SW.cpp (dlopen wrapper) is compiled into eonclib via eonclib_sources.
