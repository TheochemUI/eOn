# Fortran-only shared module: loaded at runtime via dlopen/LoadLibrary.
# Neighbour list built via the vendored vesin Fortran interface.
eon_aluminum = shared_module(
    'eon_aluminum',
           'vesin_al_helper.f90',
           'fofrhoDblexp.f',
           'dfrhoDblexp.f',
           'potinit.f',
           'forces.f',
           'gagafeDblexp.f',
           'embedenergy.f',
           'sumembforce.f',
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
                'eon_aluminum.def',
            ),
        ] : []
    ),
    install: true,
)

# Aluminum.cpp (dlopen wrapper) is compiled into eonclib via eonclib_sources.
