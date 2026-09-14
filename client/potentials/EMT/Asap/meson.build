# Headers live under include/eon/potentials/EMT/Asap/; short includes
# (Atoms.h, EMT.h, ...) need that directory on the include path.
# Path must be relative to this meson file (Meson rejects absolute
# include_directories into the source tree).
_asap_hdr_inc = include_directories(
    '../../../../include/eon/potentials/EMT/Asap',
)
asap = library(
    'asap',
    [
        'NeighborList.cpp',
        'GhostAtoms.cpp',
        'CellList.cpp',
        'Vec.cpp',
        'Atoms.cpp',
        'Exception.cpp',
        'EMTRasmussenParameterProvider.cpp',
        'EMTDefaultParameterProvider.cpp',
        'EMT.cpp',
        'stacktrace.c',
        'SuperCell.cpp',
    ],
    dependencies: _deps,
    link_with: _linkto,
    cpp_args: _args,
    include_directories: _incdirs + [_asap_hdr_inc],
    install: true,
)
