# Fortran-only shared library: loaded at runtime via dlopen
eon_sw = shared_library('eon_sw',
           'SWFortran.f90',
           fortran_args : _fargs,
           install : true)

# C++ wrapper (compiled into client, loads Fortran lib at runtime)
sw_cpp = static_library('sw_cpp',
           'SW.cpp',
           dependencies : _deps,
           cpp_args : _args,
           link_with : _linkto,
           include_directories: _incdirs)
