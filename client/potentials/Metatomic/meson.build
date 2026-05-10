mta_libdir = meson.current_source_dir()
# Precompile the torch + metatensor + metatomic header chain. Single
# heaviest TU in the project; PCH cuts cold rebuild from minutes to
# seconds. See pch/metatomic_pch.h for the rationale.
metatomic_pot = library(
    'metatomic_pot',
    'MetatomicPotential.cpp',
    cpp_pch: 'pch/metatomic_pch.h',
    dependencies: _deps,
    cpp_args: _args,
    link_with: _linkto,
    include_directories: _incdirs,
    install: true,
)
