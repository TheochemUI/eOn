parsers = library(
    'parsers',
    sources: files('ParseJob.cc', 'ParsePot.cc'),
    include_directories: _incdirs,
    dependencies: _deps,
    link_with: _linkto,
    cpp_args: _args,
    install: true,
)
_linkto += parsers
