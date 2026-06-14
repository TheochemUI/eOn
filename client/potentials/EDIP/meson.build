# Fortran-only shared library: loaded at runtime via dlopen
eon_edip = shared_library('eon_edip',
           'edipFortran.f90',
           fortran_args : _fargs,
           install : true)

# C++ wrapper (compiled into client, loads Fortran lib at runtime)
edip_cpp = static_library('edip_cpp',
           'EDIP.cpp',
           dependencies : _deps,
           cpp_args : _args,
           link_with : _linkto,
           include_directories: _incdirs)
