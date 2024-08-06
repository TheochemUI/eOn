ecmat_srcs = [
    '../Element.cc',
    'Matter.cpp',
    'MatterHelpers.cc',
    'StructComparer.cc',
]

_linkto += library(
    'ecmat',
    sources: files(ecmat_srcs),
    include_directories: _incdirs,
    dependencies: _deps,
    link_with: _linkto,
    cpp_args: _args,
    install: true,
)
