# Linked into eonclib. library() follows the project default_library, and
# a Windows DLL with no dllexport has no import library. static_library()
# stays an archive. A shared-default wheel cannot install that archive
# (libdir_static); the Windows test job already asks for static.
if is_windows
    extpot = static_library(
        'extpot',
        'ExtPot.cpp',
        dependencies: _deps,
        cpp_args: _args,
        link_with: _linkto,
        include_directories: _incdirs,
        install: get_option('default_library') != 'shared',
    )
else
    extpot = library(
        'extpot',
        'ExtPot.cpp',
        dependencies: _deps,
        cpp_args: _args,
        link_with: _linkto,
        include_directories: _incdirs,
        install: true,
    )
endif
