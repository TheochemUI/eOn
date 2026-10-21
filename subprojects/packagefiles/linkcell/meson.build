project(
    'linkcell',
    'c',
    version: '0.3.11',
    meson_version: '>=1.2.0',
    license: 'MIT',
)

# Headers only. The C ABI is linked from client/rust/eon-rust-abi.

install_headers(
    'include/linkcell.h',
    'include/linkcell.hpp',
    'include/linkcell_gpu.h',
    'include/linkcell_gpu.hpp',
)

linkcell_dep = declare_dependency(
    include_directories: include_directories('include'),
    version: meson.project_version(),
)
# The wrap provides this name. The umbrella build does not compile gpulite.
linkcell_gpu_dep = disabler()
meson.override_dependency('linkcell', linkcell_dep)
