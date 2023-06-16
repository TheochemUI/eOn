subdir('Asap')

emt = shared_library('emt',
           'EffectiveMediumTheory.cpp',
           dependencies : eon_deps,
           cpp_args : eon_extra_args,
           link_with : asap,
           install : true)
