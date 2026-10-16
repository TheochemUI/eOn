project(
    'readcon-db',
    'c',
    version: '0.1.6',
    meson_version: '>=1.3.0',
    default_options: [
        'buildtype=debugoptimized',
        'warning_level=3',
        'c_std=c99',
        'b_ndebug=if-release',
    ],
    license: 'MIT',
)

# include/readcon-db.h is hand-written and shipped. cbindgen is not used.
fs = import('fs')
_msproot = meson.project_source_root()
_mbproot = meson.project_build_root()
_incdir = include_directories('include')

if not fs.is_file(_msproot / 'include' / 'readcon-db.h')
    error(
        'missing shipped header include/readcon-db.h',
    )
endif
if not fs.is_file(_msproot / 'include' / 'readcon-db-mpi.h')
    error(
        'missing shipped header include/readcon-db-mpi.h',
    )
endif

cc = meson.get_compiler('c')
cargo_prog = find_program('cargo', required: true)
rustc_prog = find_program('rustc', required: true)
_py = find_program('python3', 'python')
_rv = run_command(rustc_prog, '--version', check: true).stdout().strip().split()
if _rv.length() < 2 or not _rv[1].version_compare('>= 1.75.0')
    error('readcon-db requires rustc >= 1.75.0')
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
    _cargo_shared = 'readcon_db.dll'
    _cargo_static = 'readcon_db.lib'
else
    _so = host_machine.system() == 'darwin' ? '.dylib' : '.so'
    _cargo_shared = 'libreadcon_db' + _so
    _cargo_static = 'libreadcon_db.a'
endif

_profile = get_option('buildtype') == 'debug' ? 'debug' : 'release'

_src_rs = run_command(
    _py,
    '-c',
    'from pathlib import Path; print("\\n".join(str(p) for p in sorted(Path("src").rglob("*.rs"))))',
    check: true,
).stdout().strip().split('\n')

# MSVC links a shared library through readcon_db.dll.lib. rustc writes
# that import library beside the DLL (often only under deps/). The
# static archive is readcon_db.lib and is a different file. A path
# string is not a ninja edge, so the import library is its own output
# and readcon_dep link_with's it.
_cargo_outputs = [_cargo_shared, _cargo_static]
_cargo_install_dirs = [get_option('libdir'), get_option('libdir')]
_cargo_cmd = [
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
]
if host_machine.system() == 'windows'
    _cargo_outputs += ['readcon_db.dll.lib']
    _cargo_install_dirs += [get_option('libdir')]
    _cargo_cmd += ['@OUTPUT2@']
endif
_cargo_cmd += [get_option('cargo_features')]

readcon_cargo = custom_target(
    'readcon-cargo',
    input: files('Cargo.toml', 'src/lib.rs'),
    depend_files: files('Cargo.lock') + files(_src_rs),
    output: _cargo_outputs,
    command: _cargo_cmd,
    build_by_default: true,
    install: true,
    install_dir: _cargo_install_dirs,
)

install_headers(
    'include/readcon-db.h',
    'include/readcon-db-mpi.h',
)

# link_args + sources: Meson must not promote an installed static consumer
# to link_whole (custom_target is not a library target). On Windows the
# shared library is linked with link_with of the import library instead,
# which is the order edge link.exe can open.
_link_idx = get_option('default_library') == 'static' ? 1 : 0
if host_machine.system() == 'windows' and get_option('default_library') != 'static'
    readcon_implib = readcon_cargo[2]
    readcon_dep = declare_dependency(
        link_with: readcon_implib,
        dependencies: _rust_sys_deps,
        include_directories: _incdir,
        version: meson.project_version(),
    )
elif host_machine.system() == 'windows'
    readcon_dep = declare_dependency(
        link_with: readcon_cargo[1],
        dependencies: _rust_sys_deps,
        include_directories: _incdir,
        version: meson.project_version(),
    )
else
    readcon_dep = declare_dependency(
        link_args: [readcon_cargo[_link_idx].full_path()],
        sources: readcon_cargo,
        dependencies: _rust_sys_deps,
        include_directories: _incdir,
        version: meson.project_version(),
    )
endif
meson.override_dependency('readcon-db', readcon_dep)

if get_option('with_tests')
    test(
        'cargo-tests',
        cargo_prog,
        args: ['test', '--manifest-path', _msproot / 'Cargo.toml', '--locked'],
        timeout: 300,
    )
endif


pkg = import('pkgconfig')
pkg.generate(
    name: 'readcon-db',
    filebase: 'readcon-db',
    description: 'Mmap-backed CON frame corpus (LMDB) with C/C++ FFI',
    version: meson.project_version(),
    # Raw `-lfoo` strings do not get `-L${libdir}`; wrap consumers then
    # fail at link with "cannot find -lreadcon_db".
    libraries: ['-L${libdir}', '-lreadcon_db'],
    url: 'https://github.com/lode-org/readcon-db',
)
