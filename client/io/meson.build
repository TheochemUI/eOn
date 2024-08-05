_linkto += library(
    'ecio_helpers',
    'IOHelpers.cc',
    include_directories: _incdirs,
    dependencies: _deps,
    link_with: _linkto,
    cpp_args: _args,
    install: true,
)

subdir('writers')

# Creator function

_linkto += library(
    'ecio_creat',
    'WriteCreator.cc',
    include_directories: _incdirs,
    dependencies: _deps,
    link_with: _linkto,
    cpp_args: _args,
    install: true,
)
