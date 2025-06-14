mta_libdir = meson.current_source_dir()
metatomic_pot = shared_library(
    'metatomic_pot',
    'MetatomicPotential.cpp',
    dependencies: _deps,
    cpp_args: _args,
    link_with: _linkto,
    include_directories: _incdirs,
    install: true,
)
