# Fortran-only shared module: loaded at runtime via dlopen/LoadLibrary.
# Neighbour list built via the vendored vesin Fortran interface.
eon_lenosky = shared_module(
    'eon_lenosky',
           'lenoskyFortran.f90',
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
                'eon_lenosky.def',
            ),
        ] : []
    ),
    install: true,
)

# Lenosky.cpp (dlopen wrapper) is compiled into eonclib via eonclib_sources.
