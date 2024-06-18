subdir('Asap')

emt = shared_library('emt',
           'EffectiveMediumTheory.cpp',
           dependencies : _deps,
           cpp_args : _args,
           link_with : [_linkto, asap],
           include_directories: _incdirs,
           install : true)
