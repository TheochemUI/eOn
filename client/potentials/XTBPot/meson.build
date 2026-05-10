# XTB potential: loaded at runtime via dlopen, no compile-time dependency.
# XtbLoader mirrors LammpsLoader: opens libxtb on first instance() call,
# caches every xtb_* function pointer, throws via require_loaded() when
# the user actually constructs an XTBPot. Skipping the link-time
# dependency('xtb') means the same eonclient binary works with or
# without xtb installed in the conda env.
xtb_eon = library('xtbeon',
           'XTBPot.cpp',
           'XtbLoader.cpp',
           dependencies : _deps,
           cpp_args : _args,
           link_with : _linkto,
           include_directories: _incdirs,
           install : true)
