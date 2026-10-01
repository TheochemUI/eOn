_rgpot_deps = _deps
_rgpot_cpp = _args
# EON_RGPOT_MPI compiles the calculator-group construction agreement and the
# grouped exit handler. It follows rgpot's MPI (-Drgpot:with_mpi=enabled, or
# an installed rgpot built with MPI), which is how calculator groups are
# built. eOn's with_mpi is the AKMC client/server switch; it still turns the
# macro on when MPI is found. A system mpic++ adds -I/usr/include, which
# breaks a sysroot compiler, so neither path probes MPI on its own.
if _rgpot_mpi_dep.found()
    _rgpot_deps += [_rgpot_mpi_dep]
    _rgpot_cpp += ['-DEON_RGPOT_MPI']
elif get_option('with_mpi').enabled()
    _rgpot_mpi = dependency('mpi', language: 'cpp', required: false)
    if _rgpot_mpi.found()
        _rgpot_deps += [_rgpot_mpi]
        _rgpot_cpp += ['-DEON_RGPOT_MPI']
    endif
endif

rgpot_pot = library(
    'rgpot_pot',
    sources: [
        'RgpotPot.cpp',
        'RGPotEngine.cpp',
        'MetatomicEngineLoader.cpp',
        'XTBEngineLoader.cpp',
    ],
    include_directories: _incdirs,
    dependencies: _rgpot_deps,
    cpp_args: _rgpot_cpp,
    link_with: _linkto,
    install: true,
)
