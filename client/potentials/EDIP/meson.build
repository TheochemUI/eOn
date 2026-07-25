# Fortran-only shared module: loaded at runtime via dlopen/LoadLibrary.
# Neighbour list built via the vendored vesin Fortran interface.
eon_edip = shared_module(
    'eon_edip',
           'edipFortran.f90',
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
                'eon_edip.def',
            ),
        ] : []
    ),
    install: true,
)

# EDIP.cpp (dlopen wrapper) is compiled into eonclib via eonclib_sources.
