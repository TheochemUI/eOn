add_languages('cython', required: true)
# Kanged from scipy
cy = meson.get_compiler('cython')
libqsc = static_library(
    'qsc',
    '_qsc.pyx',
    override_options : ['cython_language=cpp'],
    dependencies : _deps,
    include_directories: _incdirs,
    link_with : _linkto,
)
qsc=shared_library('qsc',
           'QSC.cpp',
           dependencies : _deps,
           cpp_args : _args,
           link_with : _linkto + libqsc,
           include_directories: _incdirs,
           install : true)
