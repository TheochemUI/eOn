# Fortran-only shared library: loaded at runtime via dlopen
eon_fehe = shared_library('eon_fehe',
           'feforce.f',
           fortran_args : _fargs,
           install : true)

# C++ wrapper (compiled into client, loads Fortran lib at runtime)
fehe_cpp = static_library('fehe_cpp',
           'FeHe.cpp',
           dependencies : _deps,
           cpp_args : _args,
           link_with : _linkto,
           include_directories: _incdirs)
