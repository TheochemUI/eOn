# Fortran-only shared library: loaded at runtime via dlopen
eon_tersoff = shared_library('eon_tersoff',
           'tersoffFortran.f90',
           fortran_args : _fargs,
           install : true)

# C++ wrapper (compiled into client, loads Fortran lib at runtime)
tersoff_cpp = static_library('tersoff_cpp',
           'Tersoff.cpp',
           dependencies : _deps,
           cpp_args : _args,
           link_with : _linkto,
           include_directories: _incdirs)
