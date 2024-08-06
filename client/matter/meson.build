ecmat_srcs = [
    '../Parameters.cpp',  # TODO(rg) remove
    '../Element.cc',
    'Matter.cpp',
    'MatterHelpers.cc',
    '../HelperFunctions.cpp',
]

ecmat = library(
    'ecmat',
    sources: files(ecmat_srcs),
    include_directories: _incdirs,
    dependencies: _deps,
    link_with: _linkto,
    cpp_args: _args,
    install: true,
)

_linkto += ecmat
