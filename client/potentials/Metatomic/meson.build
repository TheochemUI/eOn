mta_libdir = meson.current_source_dir()
mta_so = cppc.find_library('metatomic', dirs : mta_libdir)
metatomic_pot = shared_library('metatomic_pot',
           'MetatomicPotential.cpp',
           dependencies : [ _deps, mta_so],
           cpp_args : _args,
           link_with : _linkto,
           include_directories: _incdirs,
           install : true)
