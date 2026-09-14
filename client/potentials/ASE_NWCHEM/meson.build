# fpe_handler lives in eoncbase (_linkto); EnvHelpers is client-root source.
asenwchem = library(
    'asenwchem',
    [
        'ASE_NWCHEM.cpp',
        '../../EnvHelpers.cc',
    ],
    dependencies: _deps,
    cpp_args: _args,
    link_with: _linkto,
    include_directories: _incdirs,
    install: true,
)
