project(
    'readcon-core',
    'c',
    version: '0.14.5',
    meson_version: '>=1.3.0',
    default_options: [
        'buildtype=debugoptimized',
        'warning_level=3',
        'c_std=c99',
        'b_ndebug=if-release',
    ],
    license: 'MIT',
)

# Shipped C/C++ headers. cbindgen is maintainer-only
# (scripts/regen-capi-headers.sh); wrap and FetchContent consumers never run it.
fs = import('fs')
_msproot = meson.project_source_root()
_mbproot = meson.project_build_root()
_incdir = include_directories('include')

if not fs.is_file(_msproot / 'include' / 'readcon-core.h')
    error(
        'missing shipped header include/readcon-core.h; maintainers run scripts/regen-capi-headers.sh',
    )
endif

cc = meson.get_compiler('c')
cargo_prog = find_program('cargo', required: true)
rustc_prog = find_program('rustc', required: true)
_py = find_program('python3', 'python')
_rv = run_command(rustc_prog, '--version', check: true).stdout().strip().split()
if _rv.length() < 2 or not _rv[1].version_compare('>= 1.88.0')
    error('readcon-core requires rustc >= 1.88.0')
endif

_rust_sys_deps = []
if host_machine.system() == 'linux'
    _rust_sys_deps += [
        dependency('threads'),
        cc.find_library('dl', required: true),
        cc.find_library('m', required: true),
    ]
elif host_machine.system() == 'darwin'
    _rust_sys_deps += [dependency('threads')]
elif host_machine.system() == 'windows'
    _rust_sys_deps += [
        cc.find_library('ws2_32', required: true),
        cc.find_library('userenv', required: true),
        cc.find_library('bcrypt', required: true),
        cc.find_library('ntdll', required: true),
    ]
endif

if host_machine.system() == 'windows'
    _cargo_shared = 'readcon_core.dll'
    _cargo_static = 'readcon_core.lib'
else
    _so = host_machine.system() == 'darwin' ? '.dylib' : '.so'
    _cargo_shared = 'libreadcon_core' + _so
    _cargo_static = 'libreadcon_core.a'
endif

_profile = get_option('buildtype') == 'debug' ? 'debug' : 'release'

_src_rs = run_command(
    _py,
    '-c',
    'from pathlib import Path; print("\\n".join(str(p) for p in sorted(Path("src").rglob("*.rs"))))',
    check: true,
).stdout().strip().split('\n')

# Cargo still emits both crate types (cdylib + staticlib). Do not wrap
# the archive in a Meson static_library: that is the link_whole trap
# the 0.14 custom_target + link_args pattern exists to avoid.
readcon_cargo = custom_target(
    'readcon-cargo',
    input: files('Cargo.toml', 'src/lib.rs'),
    depend_files: files('Cargo.lock', 'build.rs') + files(_src_rs),
    output: [_cargo_shared, _cargo_static],
    command: [
        _py,
        files('scripts/meson_cargo_build.py'),
        cargo_prog.full_path(),
        _msproot,
        _mbproot / 'cargo-target',
        _profile,
        _cargo_shared,
        _cargo_static,
        '@OUTPUT0@',
        '@OUTPUT1@',
        get_option('cargo_features'),
    ],
    build_by_default: true,
    install: true,
    install_dir: [get_option('libdir'), get_option('libdir')],
)

# Indexing readcon_cargo in declare_dependency.sources does not create a
# ninja edge. Putting the cdylib in sources records NEEDED / LNK1107.
# A generated .c that depends on the whole cargo target is the order
# edge; consumers add it as a source even when they partial_dependency
# the rest of readcon_dep.
readcon_order = custom_target(
    'readcon-order',
    output: 'readcon_order.c',
    command: [
        _py,
        files('scripts/write_order_stamp.py'),
        '@OUTPUT@',
    ],
    depends: [readcon_cargo],
    build_by_default: true,
)

install_headers(
    'include/readcon-core.h',
    'include/readcon-core.hpp',
    'include/readcon-metatensor.h',
)

_link_idx = get_option('default_library') == 'static' ? 1 : 0
readcon_dep = declare_dependency(
    link_args: [readcon_cargo[_link_idx].full_path()],
    sources: readcon_order,
    dependencies: _rust_sys_deps,
    include_directories: _incdir,
    version: meson.project_version(),
)
meson.override_dependency('readcon-core', readcon_dep)

if get_option('with_tests')
    test(
        'cargo-tests',
        cargo_prog,
        args: ['test', '--manifest-path', _msproot / 'Cargo.toml', '--locked'],
        timeout: 300,
    )
endif

if get_option('with_examples')
    subdir('examples')
endif

pkg = import('pkgconfig')
pkg.generate(
    name: 'readcon-core',
    filebase: 'readcon-core',
    description: 'CON/convel file reader and writer with FFI, Python, Julia bindings',
    version: meson.project_version(),
    # Raw `-lfoo` strings do not get `-L${libdir}`; wrap consumers then
    # fail at link with "cannot find -lreadcon_core".
    libraries: ['-L${libdir}', '-lreadcon_core'],
    url: 'https://github.com/lode-org/readcon-core',
)
