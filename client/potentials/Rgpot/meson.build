_rgpot_deps = _deps
_rgpot_cpp = _args
# A system mpic++ adds -I/usr/include. That path breaks a sysroot compiler,
# so MPI is the same switch as the MPI potential.
if get_option('with_mpi')
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
