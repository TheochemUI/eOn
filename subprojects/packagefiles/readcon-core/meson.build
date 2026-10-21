project(
    'readcon-core',
    'c',
    version: '0.16.0',
    meson_version: '>=1.3.0',
    license: 'MIT',
)

# Shipped C/C++ headers. The staticlib is not built here: eOn links one
# umbrella archive (client/rust/eon-rust-abi) that depends on this crate
# as an rlib. cbindgen stays maintainer-only.

fs = import('fs')
if not fs.is_file(meson.project_source_root() / 'include' / 'readcon-core.h')
    error(
        'missing shipped header include/readcon-core.h; maintainers run scripts/regen-capi-headers.sh',
    )
endif

install_headers(
    'include/readcon-core.h',
    'include/readcon-core.hpp',
    'include/readcon-metatensor.h',
)

readcon_dep = declare_dependency(
    include_directories: include_directories('include'),
    version: meson.project_version(),
)
meson.override_dependency('readcon-core', readcon_dep)
