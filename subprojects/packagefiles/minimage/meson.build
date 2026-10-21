project(
    'minimage',
    'c',
    version: '0.1.4',
    meson_version: '>=1.2.0',
    license: 'MIT',
)

# Headers only. eOn links the C ABI from client/rust/eon-rust-abi so
# readcon-core, linkcell, and minimage share one Rust standard library.
#
# Upstream lode-org/minimage meson.build (v0.1.4) always sets
#   CARGO_TARGET_<TRIPLE>_LINKER = cc.cmd_array()[-1]
# On MSVC that value is cl. rustc links its build script with cl instead of
# link.exe and the link fails: LNK2001 memcpy, mainCRTStartup,
# __CxxFrameHandler3. The upstream fix is to set that environment variable
# only when cc.get_argument_syntax() != 'msvc'.

install_headers(
    'include/minimage.h',
    'include/minimage.hpp',
)

minimage_dep = declare_dependency(
    include_directories: include_directories('include'),
    version: meson.project_version(),
)
meson.override_dependency('minimage', minimage_dep)
