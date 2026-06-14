# Fortran-only shared library: loaded at runtime via dlopen
eon_aluminum = shared_library('eon_aluminum',
           'fofrhoDblexp.f',
           'dfrhoDblexp.f',
           'potinit.f',
           'forces.f',
           'gagafeDblexp.f',
           'embedenergy.f',
           'sumembforce.f',
           fortran_args : _fargs,
           install : true)

# C++ wrapper (compiled into client, loads Fortran lib at runtime)
aluminum_cpp = static_library('aluminum_cpp',
           'Aluminum.cpp',
           dependencies : _deps,
           cpp_args : _args,
           link_with : _linkto,
           include_directories: _incdirs)
