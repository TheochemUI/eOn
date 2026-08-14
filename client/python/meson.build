# pyeonclient — nanobind-only bindings to libeonclib
#
# Policy (no pybind11 for this module):
#   * Nanobind only for server-facing client bindings.
#   * Stable ABI (abi3 / Py_LIMITED_API 3.12) on GIL builds with Python >= 3.12.
#   * Free-threaded (NB_FREE_THREADED) when the build interpreter has
#     Py_GIL_DISABLED (CPython >= 3.13t). No stable ABI for free-threaded yet
#     (nanobind ignores STABLE_ABI in that case; we detect and choose one path).
#
# Note: pybind11 remains only for embedding ASE pots inside eonclient C++
# (with_ase*), which is a separate surface from pyeonclient.
#
# Layout:
#   pyeonclient/__init__.py
#   pyeonclient/_core*.so | _core*.abi3.so
#   pyeonclient/bridge.py

_nb_probe = run_command(
    py,
    '-c',
    '''
import nanobind, pathlib, sys, sysconfig
p = pathlib.Path(nanobind.__file__).resolve().parent
cands = [p / "src" / "nb_combined.cpp", p / "nb_combined.cpp"]
comb = next(c for c in cands if c.is_file())
# nanobind 2.x vendors tsl/robin_map under ext/robin_map/include
robin = p / "ext" / "robin_map" / "include"
if not (robin / "tsl" / "robin_map.h").is_file():
    # older layouts / system installs
    for alt in (p / "include", pathlib.Path(nanobind.include_dir())):
        if (alt / "tsl" / "robin_map.h").is_file():
            robin = alt
            break
ft = bool(sysconfig.get_config_var("Py_GIL_DISABLED"))
ver = sys.version_info
print(nanobind.include_dir())
print(robin)
print(comb)
print("1" if ft else "0")
print(f"{ver.major}.{ver.minor}")
print(nanobind.__version__)
''',
    check: false,
)
if _nb_probe.returncode() != 0
    error(
        'with_pyeonclient=true requires nanobind>=2.2 in the build Python '
        + '(pixi/pip install nanobind)',
    )
endif
_nb_lines = _nb_probe.stdout().strip().split('\n')
_nb_include = _nb_lines[0]
_nb_robin = _nb_lines[1]
_nb_combined = _nb_lines[2]
_nb_free_threaded = _nb_lines[3] == '1'
_nb_version = _nb_lines[5]

message('pyeonclient: nanobind @0@'.format(_nb_version))
message('pyeonclient: free-threaded interpreter = @0@'.format(_nb_free_threaded))
message('pyeonclient: robin_map include = @0@'.format(_nb_robin))

# Use -isystem compile_args (not include_directories) so absolute paths under a
# project-local .pixi env are allowed; Meson forbids abs include_directories
# into the source tree.
nanobind_dep = declare_dependency(
    compile_args: [
        '-isystem' + _nb_include,
        '-isystem' + _nb_robin,
    ],
)

pyeonclient_sources = files(
    'bind/module.cpp',
    'bind/bind_enums.cpp',
    'bind/bind_parameters.cpp',
    'bind/bind_potential.cpp',
    'bind/bind_matter.cpp',
    'bind/bind_ase.cpp',
    'bind/bind_jobs.cpp',
    'bind/bind_neb.cpp',
    'bind/bind_eigenmode.cpp',
    'bind/bind_saddle.cpp',
    'bind/bind_analysis.cpp',
    'bind/bind_sampling.cpp',
)

_nb_incs = [
    _incdirs,
    include_directories('..'),
    include_directories('bind'),
]

_nb_cpp_args = _args + [
    '-DNB_DOMAIN=pyeonclient',
    '-Wno-unused-parameter',
]

if _nb_free_threaded
    # Free-threaded nanobind ABI (no abi3 until PEP 803)
    _nb_cpp_args += ['-DNB_FREE_THREADED']
    message('pyeonclient: FREE_THREADED build (NB_FREE_THREADED); stable ABI skipped')
    py.extension_module(
        '_core',
        sources: pyeonclient_sources + [_nb_combined],
        include_directories: _nb_incs,
        dependencies: [_exe_deps, nanobind_dep],
        link_with: _linkto,
        cpp_args: _nb_cpp_args,
        link_args: _linkargs,
        install: true,
        subdir: 'pyeonclient',
        build_rpath: ':'.join(_runtime_rpath),
        install_rpath: ':'.join(_runtime_rpath),
        gnu_symbol_visibility: 'hidden',
    )
else
    # Stable ABI abi3 when Python >= 3.12
    _use_abi3 = py.language_version().version_compare('>=3.12')
    if _use_abi3
        # limited_api sets -DPy_LIMITED_API; do not also define it in cpp_args
        message('pyeonclient: STABLE_ABI build (abi3, Py_LIMITED_API=3.12)')
        py.extension_module(
            '_core',
            sources: pyeonclient_sources + [_nb_combined],
            include_directories: _nb_incs,
            dependencies: [_exe_deps, nanobind_dep],
            link_with: _linkto,
            cpp_args: _nb_cpp_args,
            link_args: _linkargs,
            install: true,
            subdir: 'pyeonclient',
            build_rpath: ':'.join(_runtime_rpath),
            install_rpath: ':'.join(_runtime_rpath),
            gnu_symbol_visibility: 'hidden',
            limited_api: '3.12',
        )
    else
        warning(
            'pyeonclient: Python < 3.12 — building version-specific extension '
            + '(stable ABI requires 3.12+)',
        )
        py.extension_module(
            '_core',
            sources: pyeonclient_sources + [_nb_combined],
            include_directories: _nb_incs,
            dependencies: [_exe_deps, nanobind_dep],
            link_with: _linkto,
            cpp_args: _nb_cpp_args,
            link_args: _linkargs,
            install: true,
            subdir: 'pyeonclient',
            build_rpath: ':'.join(_runtime_rpath),
            install_rpath: ':'.join(_runtime_rpath),
            gnu_symbol_visibility: 'hidden',
        )
    endif
endif

py.install_sources(
    [
        'pyeonclient/__init__.py',
        'pyeonclient/api.py',
        'pyeonclient/backends.py',
        'pyeonclient/bridge.py',
        'pyeonclient/ase_bridge.py',
        'pyeonclient/gpr_overlay.py',
        'pyeonclient/models.py',
        'pyeonclient/steps.py',
    ],
    pure: false,
    subdir: 'pyeonclient',
)
