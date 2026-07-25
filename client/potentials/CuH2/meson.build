# Neighbour list built via the vendored vesin Fortran interface.
cuh2 = library('cuh2',
           'CuH2.cpp',
           'eam_isoc.f90',
           'vesin_cuh2.f90',
           'eamroutines.f90',
           dependencies : _deps + [vesin_f_dep],
           link_with: _linkto,
           include_directories: _incdirs,
           cpp_args : _args,
           fortran_args : _fargs + ['-fno-implicit-none'],
           install : true)
