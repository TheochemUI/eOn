_rgpot_deps = _deps
_rgpot_cpp = _args
# librgpot_pot carries no MPI dependency, so a plain eonclient never maps
# libmpi and its transports. The calculator-group construction agreement
# and the grouped exit handler are librgpot_pot_mpi, which RGPotEngine
# dlopens only under an MPI launcher. It builds when MPI is available:
# rgpot's MPI (-Drgpot:with_mpi=enabled, or an installed rgpot built with
# MPI), or eOn's with_mpi. A system mpic++ adds -I/usr/include, which breaks
# a sysroot compiler, so neither path probes MPI on its own.
_rgpot_group_mpi = _rgpot_mpi_dep
if not _rgpot_group_mpi.found() and get_option('with_mpi').enabled()
    _rgpot_group_mpi = dependency('mpi', language: 'cpp', required: false)
endif

rgpot_pot = library(
    'rgpot_pot',
    sources: [
        'RgpotPot.cpp',
        'RGPotEngine.cpp',
        'GenericEngineLoader.cpp',
        'MetatomicEngineLoader.cpp',
        'XTBEngineLoader.cpp',
    ],
    include_directories: _incdirs,
    dependencies: _rgpot_deps,
    cpp_args: _rgpot_cpp,
    link_with: _linkto,
    install: true,
)

if _rgpot_group_mpi.found()
    # _Exit in the grouped handler skips the counter flush.
    _rgpot_mpi_cpp = _rgpot_cpp
    if get_option('b_coverage')
        _rgpot_mpi_cpp += ['-DEON_GCOV_DUMP']
    endif
    rgpot_pot_mpi = shared_library(
        'rgpot_pot_mpi',
        sources: ['RgpotGroupMpi.cpp'],
        include_directories: _incdirs,
        dependencies: _rgpot_deps + [_rgpot_group_mpi],
        cpp_args: _rgpot_mpi_cpp,
        link_with: _linkto,
        install: true,
    )
endif
