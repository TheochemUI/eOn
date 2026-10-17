# Linked into eonclib. A Windows DLL with no dllexport has no import
# library, so the shared build cannot resolve ExtPot.
_extpot_opts = []
_extpot_install = true
if is_windows
    _extpot_opts = ['default_library=static']
    _extpot_install = false
endif
extpot = library('extpot',
           'ExtPot.cpp',
           dependencies : _deps,
           cpp_args : _args,
           link_with : _linkto,
           include_directories: _incdirs,
           override_options: _extpot_opts,
           install : _extpot_install)
