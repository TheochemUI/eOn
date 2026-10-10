# Public headers: $src/include/eon/ (numpy/fmt style package tree).
# Sources: client/*.cpp. Both in-tree and installed form:
#   -I$prefix/include   #include "eon/Potential.h"
# thirdparty/ stays beside sources for catch2 etc.
_incdirs += [include_directories('thirdparty')]
_runtime_rpath = []

cmake = import('cmake')
fs = import('fs')

# Argum (C++20 argument parser) - single header in thirdparty

add_languages('c', native: false, required: true)
cc = meson.get_compiler('c')
cppc = meson.get_compiler('cpp')
# PGI/nvc++ cannot link_whole a static inih into a shared eonclib.
if cppc.get_id() in ['nvidia_hpc', 'pgi'] and get_option('default_library') == 'shared'
    error('nvc++ PGI linker has no link_whole; configure with -Ddefault_library=static')
endif

# eOn compiles no Fortran of its own: the kernels live in librgpot. The
# language is still detected rather than required, because a wrap build
# compiles rgpot's kernels (rgpot asks for the language itself) and the
# Windows link against them needs the flang runtime paths found below.
# A build against an installed rgpot needs no Fortran compiler at all --
# requiring one here made every such consumer supply an unused toolchain.
#
# with_fortran / with_cuh2 stay as potential selectors: they gate the
# RgpotAdapter arms in Potential.cpp, not any compilation here.
use_fortran = false
if get_option('with_fortran') or get_option('with_cuh2')
    use_fortran = add_languages('fortran', native: false, required: false)
    if use_fortran
        fc = meson.get_compiler('fortran')
        _fargs += fc.get_supported_arguments(['-fno-implicit-none'])
    endif
endif

if get_option('use_mkl')
    mkldep = dependency('mkl-dynamic-ilp64-iomp', required: true)
    add_project_arguments('-DEIGEN_USE_MKL_ALL', language: 'cpp')
    add_project_link_arguments(
        '-lmkl_intel_lp64',
        '-lmkl_core',
        '-lmkl_sequential',
        language: 'cpp',
    )
endif

# Platform detection
host_system = host_machine.system()
is_windows = host_system == 'windows'
is_mingw = is_windows and cc.get_id() == 'gcc'
# Populated when flang_rt import libs are located (Windows MSVC + flang).
_flang_rt_dir = ''

# Add conditionals
if host_system == 'darwin'
    _args += ['-DOSX=TRUE']
    # Workaround for
    # ../MinModeSaddleSearch.cpp:173:25: error: aligned allocation function of
    # type 'void *(unsigned long, enum std::align_val_t)' is only available on macOS 10.13 or newer
    add_global_arguments(['-faligned-allocation'], language: 'cpp')
endif
if is_windows
    _args += ['-DIS_WINDOWS', '-DNOMINMAX']
    # The EAM Al Fortran potential (gagafeDblexp.f) allocates ~3.6 MB of local
    # arrays on the stack (phi/phivirst dimensioned at MAXPRS=200000).  The
    # Windows default stack is only 1 MB, causing STATUS_STACK_OVERFLOW.
    # Request 16 MB to match typical Linux ulimit defaults.
    if is_mingw
        _linkargs += ['-Wl,--stack,16777216']
    elif cppc.get_id() == 'msvc'
        _linkargs += ['/STACK:16777216']
        # conda-forge flang ships flang_rt import libs under the clang resource
        # dir (e.g. Library/lib/clang/19/lib/x86_64-pc-windows-msvc/).  Meson
        # links Fortran objects with MSVC link.exe, which does not search that
        # path by default and fails with LNK1104 on FortranRuntime.dynamic.lib
        # / flang_rt.runtime.dynamic.lib. Feedstock also puts this dir on LIB.
        # flang-new reports id 'llvm-flang' or 'flang' depending on Meson.
        if use_fortran
            if fc.get_id() in ['flang', 'llvm-flang', 'flang-new']
                _flang_rt_candidates = []
                _flang_rt_libpath = run_command(
                    fc,
                    '-print-resource-dir',
                    check: false,
                )
                if _flang_rt_libpath.returncode() == 0
                    _flang_rt_res = _flang_rt_libpath.stdout().strip()
                    # Normalize D:\... to D:/... so join_paths / fs.is_* work.
                    _flang_rt_res = _flang_rt_res.replace('\\', '/')
                    _flang_rt_candidates += [
                        join_paths(_flang_rt_res, 'lib', 'x86_64-pc-windows-msvc'),
                        join_paths(_flang_rt_res, 'lib'),
                        _flang_rt_res,
                    ]
                endif
                # Fallback: conda/pixi layout under the compiler's bindir parents.
                # win-64 flang often ships import libs in Library/lib (not only under
                # the clang resource triple dir).
                _fc_path = fc.cmd_array()[0]
                _fc_dir = fs.parent(_fc_path)
                # .../Library/bin/flang-new.exe -> .../Library
                _lib_root = fs.parent(_fc_dir)
                _flang_rt_candidates += [
                    join_paths(_lib_root, 'lib'),
                    join_paths(_lib_root, 'lib', 'clang'),
                    join_paths(_lib_root, 'bin'),
                ]
                _flang_rt_dir = ''
                foreach _cand : _flang_rt_candidates
                    if _flang_rt_dir != ''
                        break
                    endif
                    if not fs.is_dir(_cand)
                        continue
                    endif
                    # Direct match or clang/<ver>/lib/x86_64-pc-windows-msvc under cand
                    foreach _probe : [
                        _cand,
                        join_paths(_cand, 'lib', 'x86_64-pc-windows-msvc'),
                    ]
                        foreach _frt : [
                            'flang_rt.runtime.dynamic.lib',
                            'FortranRuntime.dynamic.lib',
                            'flang_rt.dynamic.lib',
                        ]
                            if fs.is_file(join_paths(_probe, _frt))
                                _flang_rt_dir = _probe
                                break
                            endif
                        endforeach
                        if _flang_rt_dir != ''
                            break
                        endif
                    endforeach
                    # Walk clang version dirs one level if cand is .../lib/clang
                    if _flang_rt_dir == '' and fs.name(_cand) == 'clang'
                        # Meson has no readdir; try common versions from recent flang.
                        foreach _ver : ['21', '20', '19', '18', '17']
                            _probe = join_paths(
                                _cand,
                                _ver,
                                'lib',
                                'x86_64-pc-windows-msvc',
                            )
                            foreach _frt : [
                                'flang_rt.runtime.dynamic.lib',
                                'FortranRuntime.dynamic.lib',
                                'flang_rt.dynamic.lib',
                            ]
                                if fs.is_file(join_paths(_probe, _frt))
                                    _flang_rt_dir = _probe
                                    break
                                endif
                            endforeach
                            if _flang_rt_dir != ''
                                break
                            endif
                        endforeach
                    endif
                endforeach
                if _flang_rt_dir != ''
                    _linkargs += ['/LIBPATH:' + _flang_rt_dir]
                    message('flang_rt LIBPATH: ' + _flang_rt_dir)
                    foreach _frt : [
                        'flang_rt.runtime.dynamic.lib',
                        'FortranRuntime.dynamic.lib',
                        'flang_rt.dynamic.lib',
                    ]
                        if fs.is_file(join_paths(_flang_rt_dir, _frt))
                            _linkargs += [join_paths(_flang_rt_dir, _frt)]
                            break
                        endif
                    endforeach
                else
                    warning(
                        'flang_rt import lib dir not found; MSVC link may LNK1104',
                    )
                endif
            endif
        endif
    endif
    # Workaround used in INIFile.cpp
    # ../11.4.0/include/c++/bits/basic_string.h:1059: std::__cxx11::basic_string<_CharT, _Traits, _Alloc>::reference std::__cxx11::basic_string<_CharT, _Traits, _Alloc>::operator[](std::__cxx11::basic_string<_CharT, _Traits, _Alloc>::size_type) [with _CharT = char; _Traits = std::char_traits<char>; _Alloc = std::allocator<char>; std::__cxx11::basic_string<_CharT, _Traits, _Alloc>::reference = char&; std::__cxx11::basic_string<_CharT, _Traits, _Alloc>::size_type = long unsigned int]: Assertion '__pos <= size()' failed.
endif

add_project_arguments(
    cppc.get_supported_arguments(
        [
            '-Wno-class-memaccess',
            '-Wno-unused-function',
            '-Wno-reorder',
            '-Wno-sign-compare',
            '-Wno-unused-but-set-parameter',
            '-Wno-unused-but-set-variable',
            '-Wno-unused-parameter',
            '-Wno-unused-result',
            '-Wno-unused-variable',
        ],
    ),
    language: 'cpp',
)

if use_fortran == true
    _fargs += fc.get_supported_arguments(
        [
            # Fortran only, these are important!
            '-Wno-compare-reals',
            '-Wno-conversion',
            '-Wno-line-truncation',
            '-Wno-surprising',
            '-Wno-unused-dummy-argument',
            '-Wno-deprecated',
            '-Wno-unused-label',
            # Common
            '-Wno-maybe-uninitialized',
            '-Wno-uninitialized',
            '-Wno-unused-but-set-parameter',
            '-Wno-unused-but-set-variable',
            '-Wno-unused-parameter',
            '-Wno-unused-variable',
        ],
    )
endif

# Safer, saner default, can be turned off in release
if get_option('buildtype') != 'release'
    _args += ['-DEON_CHECKS=TRUE']
endif

# Release optimization flags.
# Do NOT pass -flto here: compile-only LTO without matching link flags breaks
# shared pot plugins (liblammps_pot etc.) with undefined std::__format symbols
# under GCC+sccache. Use Meson's -Db_lto=true when full LTO is wanted (it
# applies LTO to both compile and link).
if get_option('buildtype') == 'release' or get_option('buildtype') == 'debugoptimized'
    _opt_args = cppc.get_supported_arguments(['-O3'])
    _args += _opt_args
    if get_option('native_arch')
        _args += cppc.get_supported_arguments(['-march=native'])
    endif
    if get_option('fast_math')
        _args += cppc.get_supported_arguments(['-ffast-math'])
        if cppc.get_id() == 'msvc'
            _args += cppc.get_supported_arguments(['/fp:fast'])
        endif
    endif
endif


# Get the version, build date, and git hash
version_script = find_program('get_version.py')
version_output = run_command(version_script, check: true).stdout().strip()
version_bits = version_output.split(',')
version = version_bits[0]
build_date = version_bits[1]
git_hash = version_bits[2]
git_hash_full = run_command(
    'git',
    '-C',
    meson.project_source_root(),
    'rev-parse',
    'HEAD',
    check: false,
).stdout().strip()
# 40 hex: substring(40) is empty and index 39 exists.
hash_is_40 = git_hash_full.substring(40) == '' and git_hash_full.substring(39) != ''
if not hash_is_40
    if fs.is_dir(meson.project_source_root() / '.git')
        error('libeon_relax_engine identity needs git rev-parse HEAD (40 hex)')
    endif
    git_hash_full = git_hash
endif

# Parse semver into components
version_parts = version.split('.')
version_major = version_parts[0].to_int()
version_minor = version_parts[1].to_int()
version_patch = version_parts[2].split('+')[0].to_int()  # Remove any +git suffix

architecture = host_machine.cpu_family()
message(host_system)
message(architecture)

# Generate features string for version.h.in
features_list = [
    ['WITH_GPRD', 'GPRD'],
    ['WITH_GP_SURROGATE', 'GP Surrogate'],
    ['WITH_CATLEARN', 'CatLearn'],
    ['WITH_VASP', 'VASP'],
    ['WITH_WATER', 'Water'],
    ['WITH_AMS', 'AMS'],
    ['WITH_XTB', 'XTB'],
    ['LAMMPS_POT', 'LAMMPS'],
    ['WITH_ASE_POT', 'ASE Potentials'],
    ['EONMPI', 'MPI'],
    ['WITH_FORTRAN', 'Fortran'],
    ['CUH2_POT', 'CuH2'],
    ['WITH_ASE_ORCA', 'ASE-ORCA'],
    ['WITH_ASE_NWCHEM', 'ASE-NWChem'],
    ['WITH_METATOMIC', 'Metatomic'],
    ['WITH_SERVE_MODE', 'Serve Mode'],
    ['WITH_RGPOT', 'RGPOT direct in-process (NWChemPot/CPMDPot dlopen)'],
    ['EMBED_PYTHON', 'Embedded Python'],
    ['WITH_XTSCI', 'xtsci-optimize backend'],
    ['WITH_RGSADDLE', 'rgsaddle band and min-mode'],
]

# features_string + version.h are assembled AFTER optional with_* blocks append
# -DWITH_* to _args (metatomic, xtb, serve, artn, …). Scanning earlier under-reports
# enabled features in eonclient --features even when those flags are linked in.

# Dependencies
if not is_windows
    # libm for Unix systems
    m_dep = cppc.find_library('m', required: false)
    _deps += m_dep
    # For building with clang
    _deps += [declare_dependency(link_args: '-lstdc++')]
endif

threads_dep = dependency('threads', required: true)
_deps += threads_dep

if get_option('with_xtsci')
    xtsci_dep = dependency(
        'xtsci-optimize',
        version: '>=0.1.0',
        method: 'pkg-config',
        required: true,
    )
    _deps += xtsci_dep
    _args += ['-DWITH_XTSCI']
    eindir_dep = dependency('eindir-core', method: 'pkg-config', required: false)
    if not eindir_dep.found()
        eindir_dep = dependency('eindir', method: 'pkg-config', required: false)
    endif
    if eindir_dep.found()
        _deps += eindir_dep
        _args += ['-DWITH_EINDIR']
    endif
endif

if get_option('with_rgsaddle')
    # 0.2.0: typed surface callback, interior-only band requests.
    rgsaddle_dep = dependency(
        'rgsaddle',
        version: '>=0.2.0',
        method: 'pkg-config',
        required: true,
    )
    _deps += rgsaddle_dep
    _args += ['-DWITH_RGSADDLE']
endif

eigen_dep = dependency('eigen3', required: false)
if not eigen_dep.found()
    eigen_dep = dependency(
        'Eigen3',
        method: 'cmake',
        modules: ['Eigen3::Eigen'],
        required: true,
    )
endif
_deps += [eigen_dep]

# conda-forge ships quill.pc; source installs often only export CMake config.
quill_dep = dependency('quill', required: false)
if not quill_dep.found()
    quill_dep = dependency(
        'quill',
        method: 'cmake',
        modules: ['quill::quill'],
        required: false,
    )
endif
if not quill_dep.found()
    quill_dep = dependency(
        'Quill',
        method: 'cmake',
        modules: ['quill::quill'],
        required: false,
    )
endif
if not quill_dep.found()
    # No installed quill: build the one subprojects/quill.wrap pins, inside
    # the version range pixi.toml allows.
    quill_dep = dependency(
        'quill',
        fallback: ['quill', 'quill_dep'],
        required: true,
    )
endif
_deps += quill_dep
_args += ['-DQUILL_DISABLE_NON_PREFIXED_MACROS']

# Always use the inih wrap as a *static* library. Shared libINIReader.so is
# not on the manylinux allowlist and breaks auditwheel repair for pyeonclient.
inireader_dep = dependency('INIReader', required: false, static: true)
if not inireader_dep.found()
    inireader_dep = dependency(
        'inireader',
        fallback: ['inih', 'INIReader_dep'],
        required: true,
        static: true,
        default_options: ['default_library=static'],
    )
endif
_deps += inireader_dep

# A /usr pkg-config result is not used and the wrap is. That result is a
# pkg-config prefix of /usr, or an includedir of /usr or /usr/include.
# A real prefix (conda, EasyBuild) is kept. dependency() caches a found
# result, so the wrap is taken with subproject(), not another
# dependency() call.
nlohmann_json_dep = dependency('nlohmann_json', required: false)
_nlohmann_prefix = ''
_nlohmann_inc = ''
if nlohmann_json_dep.found()
    _nlohmann_prefix = nlohmann_json_dep.get_variable(
        pkgconfig: 'prefix',
        default_value: '',
    )
    _nlohmann_inc = nlohmann_json_dep.get_variable(
        pkgconfig: 'includedir',
        default_value: '',
    )
endif
_nlohmann_usr = (
    _nlohmann_prefix == '/usr'
    or _nlohmann_inc == '/usr'
    or _nlohmann_inc == '/usr/include'
)
if not nlohmann_json_dep.found() or _nlohmann_usr
    if _nlohmann_usr
        message(
            'nlohmann_json: a /usr pkg-config result is not used; the wrap is',
        )
    endif
    nlohmann_json_dep = subproject(
        'nlohmann_json',
        default_options: ['default_library=static'],
    ).get_variable('nlohmann_json_dep')
endif
_deps += nlohmann_json_dep

# Highway SIMD (optional, for vectorized potentials).
# highway.wrap provides the dependency name hwy (cmake target hwy).
# --force-fallback-for=hwy and --force-fallback-for=highway select that
# wrap. dependency('libhwy') does not: pkg-config can find libhwy and
# the wrap is never configured. A found hwy whose prefix is /usr is
# rejected. cmake.subproject imports hwy_list_targets. Meson
# PGICompiler (nvc++) has no get_pie_args, so that exe dies
# ninja-write. The cmake subproject's .a archives are not readable by
# MSVC link.exe, so Windows does not fall back to the wrap.
_cxx_id = cppc.get_id()
_hwy_wrap_ok = (
    host_machine.system() != 'windows'
    and _cxx_id not in ['nvidia_hpc', 'pgi']
)
hwy_dep = dependency('hwy', allow_fallback: _hwy_wrap_ok, required: false)
_hwy_prefix = ''
_hwy_inc = ''
if hwy_dep.found()
    _hwy_prefix = hwy_dep.get_variable(pkgconfig: 'prefix', default_value: '')
    _hwy_inc = hwy_dep.get_variable(
        pkgconfig: 'includedir',
        default_value: '',
    )
endif
_hwy_usr = (
    _hwy_prefix == '/usr'
    or _hwy_inc == '/usr'
    or _hwy_inc == '/usr/include'
)
_hwy_use = hwy_dep.found() and not _hwy_usr
if not _hwy_use and _hwy_wrap_ok
    if _hwy_usr
        message(
            'Highway: a /usr pkg-config result is not used; the wrap is',
        )
    else
        message('Highway: system hwy not found, trying cmake subproject...')
    endif
    hwy_opts = cmake.subproject_options()
    hwy_opts.add_cmake_defines({
        'HWY_ENABLE_TESTS': 'OFF',
        'HWY_ENABLE_EXAMPLES': 'OFF',
        'HWY_ENABLE_INSTALL': 'OFF',
        'HWY_SYSTEM_GTEST': 'OFF',
        'BUILD_TESTING': 'OFF',
        'HWY_ENABLE_CONTRIB': 'ON',
    })
    hwy_proj = cmake.subproject('highway', required: false, options: hwy_opts)
    if hwy_proj.found()
        hwy_dep = hwy_proj.dependency('hwy')
        _hwy_use = true
    endif
endif
if _hwy_use
    _deps += hwy_dep
    _args += ['-DWITH_HIGHWAY']
    message('Highway SIMD: enabled')
else
    message('Highway SIMD: disabled (not found)')
endif

# dl_dep needed for dlopen (LammpsLoader, rgpot engine loaders)
dl_dep = cppc.find_library('dl', required: false)
_deps += dl_dep

# rgpot resolution starts here because the vesin choice below depends on
# it: a wrap-built rgpot carries its own patched vesin, and linking eOn's
# copy next to it would define every vesin symbol twice. The roles the
# resulting objects play are wired further down, once the pots exist.
#
# dependency('rgpot') with default allow_fallback configures the wrap once
# with *default* options before eOn can pass its own; get_variable() then
# fails for anything the default profile skipped. Force pkg-config-only.
_rgpot_proj = []
_have_rgpot_wrap = false
_rgpot_pkg = dependency(
    'rgpot',
    version: '>=3.4.0',
    method: 'pkg-config',
    required: false,
    allow_fallback: false,
)
if not _rgpot_pkg.found()
    # One profile must serve every eOn consumer in this configure: pots
    # always. The RPC stack is on for the direct RGPOT arm (non-Windows)
    # and for serve mode. Plain Windows builds keep capnp-rpc out.
    rgpot_opts = {
        'with_rpc_client_only': false,
        'pure_lib': true,
        'with_tests': false,
        'with_examples': false,
        # engines stay runtime dlopen. Linked XTBPot is on only when eOn
        # itself asked for -Dwith_xtb=true (the factory arm is rgpot::XTBPot).
        'with_xtb': get_option('with_xtb'),
        'with_metatomic': false,
        # ExprPot is a separate option. The basic build must compile it
        # or PotType.EXPR falls through to the generic terminating error.
        'with_expr': true,
        # SW, EDIP, Lenosky, Tersoff, EAM-Al, FeHe, CuH2, and TIP4P-H are
        # Fortran kernels inside librgpot.
        'with_fortran_pots': 'enabled',
        'with_water': get_option('with_water'),
    }
    if get_option('with_rgpot_mpi').enabled()
        rgpot_opts += {'with_mpi': 'enabled'}
    endif
    # Direct RGPOT is mandatory off Windows; serve still opts into RPC there.
    if not is_windows or get_option('with_serve')
        rgpot_opts += {'with_rpc': true}
    endif
    _rgpot_proj = subproject('rgpot', default_options: rgpot_opts)
    # Wrap-built librgpot is not installed (rgpot install: not
    # is_subproject; CI uses --skip-subprojects). eonclient's rpath
    # looks in prefix/lib for librgpot.3.dylib.
    meson.add_install_script(
        find_program(meson.project_source_root() / 'scripts' / 'install_wrap_rgpot.py'),
        _rgpot_proj.get_variable('rgpot').full_path(),
        get_option('libdir'),
    )
    _have_rgpot_wrap = true
endif

# The compiler records the prefix libdir as DT_RPATH ahead of the $ORIGIN
# entries Meson adds for each target. DT_RPATH is searched in flag order,
# and ahead of LD_LIBRARY_PATH, so a library installed into that prefix is
# the one a later build loads. Project link arguments are emitted before
# those external flags. The entries are the client build layout at every
# directory depth a client target is linked from: each potential directory,
# the client directory itself, and a wrap-built librgpot. They are $ORIGIN
# relative, so install does not keep a build directory. A path that is not
# present is an extra lookup.
if host_system == 'linux'
    _pot_out = run_command(
        find_program('find'),
        meson.current_source_dir() / 'potentials',
        '-name',
        'meson.build',
        '-printf',
        '%P\\n',
        check: true,
    )
    _pot_rels = []
    foreach _pot_line : _pot_out.stdout().splitlines()
        if _pot_line.endswith('/meson.build')
            _pot_rels += 'potentials' / _pot_line.replace('/meson.build', '')
        endif
    endforeach
    _pot_sorted = []
    foreach _pot_rel : _pot_rels
        _pot_placed = false
        _pot_next = []
        foreach _pot_old : _pot_sorted
            if not _pot_placed and _pot_rel < _pot_old
                _pot_next += _pot_rel
                _pot_placed = true
            endif
            _pot_next += _pot_old
        endforeach
        if not _pot_placed
            _pot_next += _pot_rel
        endif
        _pot_sorted = _pot_next
    endforeach
    # '.' is the client build directory (libeonclib and the executables).
    # python/ is one level under it. Nested potentials set the deeper levels.
    _origin_deps = ['.'] + _pot_sorted
    _origin_depth = 1
    foreach _pot_rel : _pot_sorted
        _pot_n = _pot_rel.split('/').length()
        if _pot_n > _origin_depth
            _origin_depth = _pot_n
        endif
    endforeach
    if _have_rgpot_wrap
        _rgpot_libdir = fs.parent(_rgpot_proj.get_variable('rgpot').full_path())
        _rgpot_rel = fs.relative_to(_rgpot_libdir, meson.current_build_dir())
        if not _rgpot_rel.startswith('/')
            _origin_deps += _rgpot_rel
        endif
    endif
    _origin_entries = []
    foreach _origin_level : range(_origin_depth + 1)
        _origin_ups = ''
        foreach _origin_step : range(_origin_level)
            _origin_ups += '../'
        endforeach
        foreach _origin_dep : _origin_deps
            if _origin_dep == '.'
                _origin_rel = _origin_ups
            else
                _origin_rel = _origin_ups + _origin_dep
            endif
            if _origin_rel.endswith('/')
                _origin_rel = _origin_rel.substring(0, -1)
            endif
            if _origin_rel == ''
                _origin_entries += '$ORIGIN'
            else
                # String join, not the / operator: path join collapses
                # '$ORIGIN/..' and drops the token.
                _origin_entries += '$ORIGIN/' + _origin_rel
            endif
        endforeach
    endforeach
    add_project_link_arguments(
        '-Wl,-rpath,' + ':'.join(_origin_entries),
        language: ['c', 'cpp'],
    )
    message(
        'build-tree rpath entries before the prefix libdir:',
        _origin_entries.length(),
    )
endif

# vesin: single neighbor-list backend for Metatomic and eonc::VesinNeighbors.
# Prefer pkg-config / system; else the rgpot subproject's vendored copy;
# else eOn's own thirdparty/vesin single TU.
# The vendored TU stays. Metatomic calls vesin_neighbors directly, and
# eonc::VesinNeighbors is a public header, so Metatomic is not the only
# consumer. An installed rgpot does not export vesin headers. The wrapper
# stays in eOn; moving it would be an rgpot API change.
vesin_dep = dependency('vesin', required: false)
if not vesin_dep.found()
    vesin_dep = dependency('vesin', method: 'cmake', required: false)
endif
if not vesin_dep.found()
    vesin_dep = cppc.find_library('vesin', required: false)
endif
if not vesin_dep.found() and is_windows
    vesin_dep = cppc.find_library('libvesin', required: false)
endif
if not vesin_dep.found() and not _rgpot_pkg.found()
    # The rgpot subproject builds vesin already, with the local patches its
    # Fortran interface binds to.
    vesin_dep = _rgpot_proj.get_variable('vesin_dep')
    message('vesin: from the rgpot subproject')
elif not vesin_dep.found()
    # vesin 0.6 threads its builds with std::thread only — do NOT add
    # OpenMP here: an unused libgomp costs ~0.3ms of load/relocation on
    # every eonclient start.
    vesin_int_deps = [dependency('threads')]
    # Static: one fewer shared object to relocate at eonclient startup
    # (the ASV point fixture is startup-dominated).
    vesin_int = static_library(
        'vesin_internal',
        ['thirdparty/vesin/vesin-single-build.cpp'],
        dependencies: vesin_int_deps,
        # Consumers link statically; hiding the symbols keeps them out of
        # the enclosing shared objects' dynamic tables (startup relocations).
        cpp_args: cppc.get_supported_arguments(['-fvisibility=hidden']),
        pic: true,
    )
    vesin_dep = declare_dependency(
        link_with: vesin_int,
        dependencies: vesin_int_deps,
        include_directories: include_directories('thirdparty/vesin'),
    )
    message('vesin: using vendored thirdparty/vesin')
else
    # Ensure #include "vesin.h" works when only a library was found.
    _incdirs += [include_directories('thirdparty/vesin')]
    message('vesin: using system / pkg-config')
endif
_deps += [vesin_dep]

eoncbase_sources = [
    'fpe_handler.cpp',
    'PotRegistry.cpp',
    'ParametersAccess.cpp',
    'PotentialParams.cpp',
    'VesinNeighbors.cpp',
    'potentials/PluginLoader.cpp',
]

# Always static. A shared eoncbase cannot see Parameters::potential_options,
# which lives in the client library, and a static install path is
# {libdir_static}, which the wheel builder cannot map.
eoncbase = static_library(
    'eoncbase',
    sources: eoncbase_sources,
    dependencies: _deps,
    cpp_args: _args,
    include_directories: _incdirs,
    pic: true,
    install: false,
)

_linkto += eoncbase

# readcon-core: Rust-based .con file I/O (required).
# Resolve before plugin subdirs. SocketNWChem includes <readcon-core.hpp>
# and links rkr_*; library() snapshots dependencies at call time, so a
# later _deps += readcon_dep does not reach libsocket_nwchem.so.
# Prefer cargo-c / pkg-config install (conda-forge host package); fall back
# to the meson subproject when the system package is not available.
# ConFileIO needs 0.16 for set_spreads_from_flat, on top of 0.15's
# displacements section (Hessian modes.con) and 0.13's clone, bulk
# set_*_from_flat and compression_from_extension.
readcon_dep = dependency(
    'readcon-core',
    version: '>=0.16.0',
    method: 'pkg-config',
    required: false,
)
readcon_order_srcs = []
# Directory of the readcon DLL this build linked. A pixi env can put
# another readcon_core.dll earlier on PATH; Windows then starts the
# test with STATUS_ENTRYPOINT_NOT_FOUND (0xC0000139).
_readcon_runtime_dir = ''
_have_readcon_dll_beside = false
if not readcon_dep.found()
    readcon_proj = subproject(
        'readcon-core',
        version: '>=0.16.0',
        default_options: {
            'with_tests': false,
            'with_examples': false,
            'default_library': 'static',
        },
    )
    readcon_dep = readcon_proj.get_variable('readcon_dep')
    # Order stamp: cargo still emits cdylib+staticlib. partial_dependency
    # drops declare_dependency.sources, so attach the stamp here. A
    # readcon-core without the stamp contributes no sources.
    readcon_order_srcs += readcon_proj.get_variable('readcon_order', [])
    if is_windows
        _readcon_dll = readcon_proj.get_variable('readcon_cargo')[0]
        _readcon_runtime_dir = fs.parent(_readcon_dll.full_path())
        # The loader searches the executable directory before PATH.
        # A pixi readcon_core.dll on PATH is a different build and the
        # process exits 0xC0000139 before main. The copy sits next to
        # test_relax_engine.exe and eon_relax_engine.dll.
        _readcon_dll_beside = custom_target(
            'readcon-dll-beside-tests',
            input: _readcon_dll,
            output: 'readcon_core.dll',
            command: [
                py,
                '-c',
                'import shutil, sys; shutil.copyfile(sys.argv[1], sys.argv[2])',
                '@INPUT@',
                '@OUTPUT@',
            ],
            build_by_default: true,
        )
        _have_readcon_dll_beside = true
    endif
endif
if is_windows
    # On Windows, do not put cargo .lib/.dll on eonclib: llvm-ar treats
    # the .lib as an archive name, and link rejects the DLL (LNK1107).
    # The order stamp is attached on eclib sources; link_args stay on
    # executables via _exe_deps. cargo-c/pkg-config has no custom_target.
    # SocketNWChem is not built on Windows.
    readcon_compile_dep = readcon_dep.partial_dependency(
        includes: true,
        compile_args: true,
    )
    _deps += [readcon_compile_dep]
else
    _deps += [readcon_dep]
endif

# readcon-db: corpus mirror. Resolve it the same way as readcon-core
# (pkg-config, then the wrap). The shared library is not archived beside
# the static readcon-core.
readcon_db_dep = dependency(
    'readcon-db',
    version: '>=0.1.6',
    method: 'pkg-config',
    required: false,
)
if not readcon_db_dep.found()
    readcon_db_proj = subproject(
        'readcon-db',
        version: '>=0.1.6',
        default_options: {
            'with_tests': false,
            'default_library': 'shared',
        },
    )
    _readcon_db_inner = readcon_db_proj.get_variable('readcon_dep')
    _readcon_db_so = readcon_db_proj.get_variable('readcon_cargo')[0]
    _readcon_db_version = _readcon_db_inner.version()
    _runtime_rpath += [fs.parent(_readcon_db_so.full_path())]
    # full_path() in link_args is a string, so ninja has no edge to the
    # cargo target and a parallel link runs while that library is absent.
    # link_with is the order edge. The --no-as-needed path stays: the
    # mirror resolves rkrdb_open at load time, and --as-needed drops a
    # library with no link-time use.
    # ld64 has no --as-needed and keeps every dylib it is given, so macOS
    # takes the bare path.
    _readcon_db_link = [_readcon_db_so.full_path()]
    if host_system != 'darwin'
        _readcon_db_link = [
            '-Wl,--no-as-needed',
            _readcon_db_so.full_path(),
            '-Wl,--as-needed',
        ]
    endif
    if is_windows
        # The mirror is compiled out on Windows (ReadconDbMirror.cpp guards
        # every rkrdb_* use with _WIN32), so nothing links readcon-db there.
        # link.exe cannot take the DLL, and without a link_with edge the
        # path is also absent when the first DLL links (LNK1181).
        readcon_db_dep = declare_dependency(
            compile_args: ['-DEON_READCON_DB_VERSION="' + _readcon_db_version + '"'],
        )
    else
        readcon_db_dep = declare_dependency(
            link_with: [_readcon_db_so],
            link_args: _readcon_db_link,
            sources: [_readcon_db_so],
            compile_args: ['-DEON_READCON_DB_VERSION="' + _readcon_db_version + '"'],
        )
    endif
else
    _readcon_db_version = readcon_db_dep.version()
    if _readcon_db_version == ''
        _readcon_db_version = '0.1.6'
    endif
    # A wrap fallback is already found(), and its link line is the cdylib
    # path. A bare -lreadcon_db has no -L under the wheel's conda linker.
    # An installed pkg-config file is kept by that same path when libdir
    # names the shared library.
    _readcon_db_keep = []
    _pc_libdir = readcon_db_dep.get_variable(
        pkgconfig: 'libdir',
        default_value: '',
    )
    if _pc_libdir != ''
        foreach _name : ['libreadcon_db.so', 'libreadcon_db.dylib']
            _cand = _pc_libdir / _name
            if fs.exists(_cand)
                _readcon_db_keep = [_cand]
                if host_system != 'darwin'
                    _readcon_db_keep = [
                        '-Wl,--no-as-needed',
                        _cand,
                        '-Wl,--as-needed',
                    ]
                endif
                break
            endif
        endforeach
    endif
    readcon_db_dep = declare_dependency(
        dependencies: [readcon_db_dep],
        link_args: _readcon_db_keep,
        compile_args: ['-DEON_READCON_DB_VERSION="' + _readcon_db_version + '"'],
    )
endif
_deps += [readcon_db_dep]

# Potentials

subdir('potentials/ExtPot')
if not is_windows
    subdir('potentials/SocketNWChem')
endif

# rgpot roles (see docs/user_guide/rgpot_integration.md):
#   non-Windows -> nwchempot_dep + cpmdpot_dep (direct dlopen engines in eOn)
#   with_serve  -> ptlrpc_dep (eOn is Cap'n Proto RPC server)
# potserv RPC *client* is not wired as an eOn Potential.
#
# Prefer an *installed* rgpot (pkg-config) for distro/EasyBuild packaging;
# fall back to subprojects/rgpot.wrap for pixi/dev. Wrap defaults (document
# for matching standalone installs — no torch/xtb at link time):
#   with_rpc_client_only=true, pure_lib=true, with_tests=false,
#   with_examples=false  (+ with_xtb=false, with_metatomic=false for packages)
ptlrpc_dep = disabler()
nwchempot_dep = disabler()
cpmdpot_dep = disabler()
_rgpot_from_pkg = false
# rgpot is a hard dependency: the classical pair pots (LJ, LJCluster,
# Morse, ZBL, ...) live in librgpot and eOn consumes them through the
# RgpotAdapter. The direct NWChem/CPMD arm is on everywhere except
# Windows. Serve mode stays optional.
#
# 1) installed module first (never wrap fallback here); _rgpot_pkg is
# probed above, next to the vesin choice that depends on it.
if _rgpot_pkg.found()
    message('rgpot: using pkg-config (installed), not subproject wrap')
    _rgpot_from_pkg = true
    # Calculator groups need rgpot's MPI, not eOn's with_mpi (the AKMC
    # client/server switch). An rgpot.pc that defines RGPOT_HAS_MPI says so.
    _rgpot_mpi_dep = dependency('', required: false)
    if cppc.get_define('RGPOT_HAS_MPI', dependencies: _rgpot_pkg) != ''
        _rgpot_mpi_dep = dependency('mpi', language: 'cpp', required: true)
    endif
    rgpot_dep = _rgpot_pkg
    nwchempot_dep = _rgpot_pkg
    cpmdpot_dep = _rgpot_pkg
    ptlrpc_dep = _rgpot_pkg
    # An installed rgpot can be built without the Fortran kernels, and then
    # SW, EDIP, Lenosky, Tersoff, EAM-Al, FeHe, CuH2, and TIP4P-H resolve to
    # nothing at link time. rgpot.pc advertises them; say so plainly here
    # rather than failing on undefined vtables.
    if cppc.get_define('RGPOT_HAS_FORTRAN_POTS', dependencies: _rgpot_pkg) == ''
        error(
            'installed rgpot was built without Fortran potentials '
            + '(rgpot.pc does not define RGPOT_HAS_FORTRAN_POTS). '
            + 'Rebuild rgpot with -Dwith_fortran_pots=enabled, or point '
            + 'PKG_CONFIG_PATH at one that has them.',
        )
    endif
    if get_option('with_water') and cppc.get_define(
        'RGPOT_HAS_WATER',
        dependencies: _rgpot_pkg,
    ) == ''
        error(
            'installed rgpot was built without the flexible water kernels '
            + '(rgpot.pc does not define RGPOT_HAS_WATER). '
            + 'Rebuild rgpot with -Dwith_water=true.',
        )
    endif
else
    # 2) wrap / multi-source fallback, configured next to the vesin choice
    # above; here it only takes on its eOn roles.
    message('rgpot: pkg-config not found; falling back to subproject wrap')
    rgpot_dep = _rgpot_proj.get_variable('rgpot_dep')
    # The subproject links MPI when -Drgpot:with_mpi=enabled; calculator
    # groups run with that and eOn's with_mpi off.
    _rgpot_mpi_dep = _rgpot_proj.get_variable(
        '_mpi_dep',
        dependency('', required: false),
    )
    # meson applies subproject default_options only on the first configure
    # of a build directory, so a directory configured before this profile
    # changed silently keeps the old values and drops potentials. Fail
    # here instead, naming the fix.
    if not _rgpot_proj.get_variable('has_fortran')
        error(
            'rgpot subproject built without Fortran potentials. '
            + 'Reconfigure with -Drgpot:with_fortran_pots=enabled, '
            + 'or wipe the build directory.',
        )
    endif
    if get_option('with_water') and not _rgpot_proj.get_variable('has_water')
        error(
            'rgpot subproject built without the flexible water kernels. '
            + 'Reconfigure with -Drgpot:with_water=true, '
            + 'or wipe the build directory.',
        )
    endif
    # pkg-config puts RGPOT_HAS_* in the .pc cflags. The wrap dependency
    # does not, so the expression potential is compiled out and
    # makePotential throws "Terminating".
    if _rgpot_proj.get_variable('has_expr')
        _args += ['-DRGPOT_HAS_EXPR']
    endif
    if get_option('with_serve')
        ptlrpc_dep = _rgpot_proj.get_variable('ptlrpc_dep')
    endif
    if not is_windows
        nwchempot_dep = _rgpot_proj.get_variable('nwchempot_dep')
        cpmdpot_dep = _rgpot_proj.get_variable('cpmdpot_dep')
    endif
endif
_deps += [rgpot_dep]
if get_option('with_serve')
    serve_capnp_dep = dependency('capnp-rpc', required: true)
    if not ptlrpc_dep.found()
        error('with_serve requires ptlrpc from rgpot (pkg-config or subproject)')
    endif
endif
if not is_windows
    if not nwchempot_dep.found() or not cpmdpot_dep.found()
        error('RGPOT requires nwchempot/cpmdpot from rgpot (pkg-config or subproject)')
    endif
    _deps += [nwchempot_dep, cpmdpot_dep]
    _args += ['-DWITH_RGPOT']
    subdir('potentials/Rgpot')
endif

# LJ, LJCluster, Morse, ZBL, and EMT live in librgpot (RgpotAdapter arms).
potentials = [
    extpot,
]
if not is_windows
    potentials += [socket_nwchem]
endif
if not is_windows
    potentials += [rgpot_pot]
endif

eonclib_sources = [
    'Parameters.cpp',
    'ParametersINI.cpp',
    'ParametersJSON.cpp',
    'ParametersSSOT.cpp',
    'ConFileIO.cpp',
    'ReadconDbMirror.cpp',
    'GeometryAnalysis.cpp',
    'Optimizer.cpp',
    'IDPPObjectiveFunction.cpp',
    'PrefactorJob.cpp',
    'LBFGS.cpp',
    'ReplicaExchangeJob.cpp',
    'BondBoost.cpp',
    'Job.cpp',
    'JobRegistry.cpp',
    'Runtime.cpp',
    'TestJob.cpp',
    'OHTSTJob.cpp',
    'GleThermostat.cpp',
    'GlobalOptimization.cpp',
    'LowestEigenmode.cpp',
    'EigenmodeStrategy.cpp',
    'MinModeSaddleSearch.cpp',
    'StructureComparisonJob.cpp',
    'SteepestDescent.cpp',
    'ImprovedDimer.cpp',
    'LORRotation.cpp',
    'PointJob.cpp',
    'Prefactor.cpp',
    'ConjugateGradients.cpp',
    'Matter.cpp',
    'PbcSimd.cpp',
    'ForceNorm.cpp',
    'FiniteDifferenceJob.cpp',
    'Lanczos.cpp',
    'Davidson.cpp',
    'MobileAtoms.cpp',
    'HessianJob.cpp',
    'InstantonJob.cpp',
    'ReplicaDynamicsJob.cpp',
    'TADJob.cpp',
    'ProcessSearchJob.cpp',
    'Bundling.cpp',
    'NEBInitialPaths.cpp',
    'NEBForceProjection.cpp',
    'NEBObjectiveFunction.cpp',
    'NEBProjection.cpp',
    'NEBOcinebController.cpp',
    'NEBSplineExtrema.cpp',
    'Tunneling.cpp',
    'QuantumFreeEnergy.cpp',
    'NEBSpringForce.cpp',
    'SolidStateNEB.cpp',
    'NEBTangent.cpp',
    'NEBZoom.cpp',
    'NudgedElasticBand.cpp',
    'MonteCarloJob.cpp',
    'DynamicsJob.cpp',
    'MonteCarlo.cpp',
    'Hessian.cpp',
    'NudgedElasticBandJob.cpp',
    'DynamicsSaddleSearch.cpp',
    'HelperFunctions.cpp',
    'RandomNumbers.cpp',
    'StringHelpers.cc',
    # MatrixHelpers.hpp is header-only under include/eon/ (via -I include).
    'Dimer.cpp',
    'Dynamics.cpp',
    'PathIntegral.cpp',
    'PIQTST.cpp',
    'PIQTSTJob.cpp',
    'GlobalOptimizationJob.cpp',
    'BiasedGradientSquaredDescent.cpp',
    'SafeHyperJob.cpp',
    'MinimizationJob.cpp',
    'JobResultProvenance.cpp',
    'Quickmin.cpp',
    'ParallelReplicaJob.cpp',
    'Potential.cpp',
    'SurrogatePotential.cpp',  # Part of the interface
    'BasinHoppingJob.cpp',
    'FIRE.cpp',
    'EpiCenters.cpp',
    'SaddleSearchJob.cpp',
    'BasinHoppingSaddleSearch.cpp',
    # Stubs when -Dwith_ira=false; NudgedElasticBand always calls IRACompare.
    'IRACompare.cpp',
    # Runtime loader always compiled; compare path stubs without WITH_IRA.
    'libs/IRA/IRAResource.cpp',
    # Runtime loader always compiled; saddle path stubs without WITH_ARTN.
    'ARTnSaddleSearch.cpp',
    'libs/ARTn/ARTnResource.cpp',
]

eonclib_sources += ['XtsciEindir.cpp']
if get_option('with_xtsci')
    eonclib_sources += ['XtsciOptimizer.cpp']
endif

if get_option('with_rgsaddle')
    eonclib_sources += ['XtsciBand.cpp', 'XtsciMinMode.cpp']
endif

eonclient_sources = ['ClientEON.cpp', 'CommandLine.cpp']

# C ABI host always compiled so Catch2 can inject IMetatomicLoader without
# the fat Torch-linked MetatomicPotential. Production still uses instance()
# when -Dwith_metatomic=true builds libmetatomic_pot.
eonclib_sources += [
    'potentials/Metatomic/MetatomicLoader.cpp',
    'potentials/Metatomic/MetatomicDynPot.cpp',
]

# ------------------- Options

# gpr_optim (+ nested rgpot) must not land in the shared pot-plugin dep set:
# those plugins link with --no-undefined and would inherit HDF5/zlib (and
# rgpot) link flags they do not need. Attach only to eonclib / executables.
#
# A checkout at subprojects/gpr_optim turns the GP dimer on under
# -Dwith_gprd=auto (the default) on Linux:
#   rsync -a ../gpr_optim/ subprojects/gpr_optim/
# Without one the build leaves min_mode_method=gprdimer off.
_gprd_opt = get_option('with_gprd')
_gprd_tree = fs.exists(meson.project_source_root() / 'subprojects' / 'gpr_optim' / 'meson.build')
if _gprd_opt.enabled() and not _gprd_tree
    error('with_gprd=enabled needs a gpr_optim checkout in subprojects/gpr_optim (rsync -a ../gpr_optim/ subprojects/gpr_optim/).')
endif
# auto: try the subproject and leave the dimer off when it cannot be
# fetched. enabled: the same failure stops configuration.
_gprd_try = _gprd_opt.enabled() or (
    _gprd_opt.auto() and host_system == 'linux' and _gprd_tree
)
_gprd_extra_deps = []
with_gprd = false
if _gprd_try
    gprd_opts = {
        'default_library': 'static',
        'use_mkl': get_option('use_mkl'),
        # AtomicGPDimer does not use gpr_optim HDF5 logs; keep HDF5 off so
        # gpr_optim_deps do not require -lz in lean CI envs.
        'use_hdf5': false,
        # The gpr_optim default looks up an installed rgmin and stops when
        # that package is absent. Disabled keeps the C++ solvers.
        'use_rgmin': 'disabled',
        'do_install': true,
        'buildtype': 'release',
        'linalg_backend': get_option('gprd_linalg_backend'),
    }
    libgprd_proj = subproject(
        'gpr_optim',
        default_options: gprd_opts,
        required: _gprd_opt,
    )
    if libgprd_proj.found()
        with_gprd = true
        libgprd = libgprd_proj.get_variable('gpr_optim_dep')
        gprd_deps = libgprd_proj.get_variable('gpr_optim_deps')
        # subdir('potentials/GPRPotential')
        # potentials += [ gprpot ]
        _args += ['-DWITH_GPRD']
        eonclib_sources += ['AtomicGPDimer.cpp', 'GPRHelpers.cpp']
        _gprd_extra_deps = [libgprd, gprd_deps]
    elif _gprd_opt.enabled()
        error('with_gprd=enabled: subproject gpr_optim did not configure.')
    else
        message('with_gprd=auto: subproject gpr_optim did not configure. GP dimer disabled.')
    endif
endif

if get_option('with_gp_surrogate')
    _args += ['-DWITH_GP_SURROGATE']
    eonclib_sources += ['GPSurrogateJob.cpp']
    if get_option('with_catlearn')
        # TODO: Cleanup, used for ase_orca too
        # Embedding the interpreter
        pyb11f_deps = [dependency('pybind11')]
        _deps += [pyb11f_deps]
        subdir('potentials/CatLearnPot')
        potentials += [catlearnpot]
        _args += ['-DWITH_CATLEARN']
    endif
endif

if get_option('with_vasp')
    if host_system != 'windows'
        subdir('potentials/VASP')
        potentials += [vasp]
        _args += ['-DWITH_VASP']
    endif
endif

if get_option('with_water')
    # The kernels live in librgpot. This flag only selects the factory arms.
    _args += ['-DWITH_WATER']
endif

if get_option('with_ams')
    # Optional, only for AMS
    absl_dep = dependency(
        'absl',
        method: 'cmake',
        modules: [
            'absl::base',
            'absl::strings',
            'absl::numeric',
            'absl::utility',
        ],
        required: true,
    )
    boost_ams_dep = dependency(
        'boost',
        modules: ['headers', 'asio', 'iostreams'],
        required: true,
    )
    ams_dep = [boost_ams_dep, absl_dep]
    _deps += ams_dep
    # End AMS
    subdir('potentials/AMS')
    subdir('potentials/AMS_IO')
    potentials += [ams, ams_io]
    _args += ['-DWITH_AMS']
endif

if get_option('with_xtb')
    # Headers for rgpot/XTBPot/XTBPot.hpp. The kernel itself is inside
    # librgpot (rgpot -Dwith_xtb=true); eOn does not compile a second one.
    # conda-forge win-64 xtb is MinGW and ships no xtb.lib. MSVC still
    # needs an import library written from the DLL.
    if host_machine.system() == 'windows' and cppc.get_id() == 'msvc'
        # conda-forge win-64 xtb is MinGW: libxtb-6.dll and libxtb.dll.a.
        # dependency('xtb') asks link.exe for xtb.lib, which that package
        # does not ship. Write the import library from the DLL.
        _xtb_pc = dependency('xtb', method: 'pkg-config', required: true)
        _xtb_prefix = _xtb_pc.get_variable(pkgconfig: 'prefix')
        _xtb_libdir = _xtb_pc.get_variable(pkgconfig: 'libdir')
        _xtb_dll = ''
        foreach _cand : [
            _xtb_prefix / 'bin' / 'libxtb-6.dll',
            _xtb_prefix / 'bin' / 'libxtb.dll',
            _xtb_prefix / 'Library' / 'bin' / 'libxtb-6.dll',
            _xtb_prefix / 'Library' / 'bin' / 'libxtb.dll',
            _xtb_libdir / 'libxtb-6.dll',
            _xtb_libdir / 'libxtb.dll',
            _xtb_prefix / 'bin' / 'xtb.dll',
        ]
            if _xtb_dll == '' and fs.exists(_cand)
                _xtb_dll = _cand
            endif
        endforeach
        if _xtb_dll == ''
            error('MSVC xtb: libxtb DLL not found under ' + _xtb_prefix)
        endif
        _xtb_machine = 'x64'
        if host_machine.cpu_family() == 'x86'
            _xtb_machine = 'X86'
        elif host_machine.cpu_family() == 'aarch64'
            _xtb_machine = 'ARM64'
        endif
        _xtb_implib = meson.current_build_dir() / 'xtb.lib'
        run_command(
            py,
            files('../scripts/ci/msvc_xtb_implib.py'),
            _xtb_dll,
            _xtb_implib,
            _xtb_machine,
            check: true,
        )
        xtb_dep = declare_dependency(
            compile_args: ['-I' + _xtb_pc.get_variable(pkgconfig: 'includedir')],
            # rgpot's xtb dependency still contributes a bare xtb.lib.
            # link.exe searches LIBPATH and the build root, not this
            # directory, so the generated import library is invisible
            # unless its directory is on the link line.
            link_args: ['/LIBPATH:' + meson.current_build_dir()],
            dependencies: cppc.find_library(
                'xtb',
                dirs: [meson.current_build_dir()],
                required: true,
            ),
            version: _xtb_pc.version(),
        )
    else
        xtb_dep = dependency('xtb', required: true)
    endif

    _deps += [xtb_dep]
    if _rgpot_pkg.found()
        if cppc.get_define('RGPOT_HAS_XTB_POT', dependencies: _rgpot_pkg) == ''
            error(
                'installed rgpot was built without XTBPot '
                + '(rgpot.pc does not define RGPOT_HAS_XTB_POT). '
                + 'Rebuild rgpot with -Dwith_xtb=true, or point '
                + 'PKG_CONFIG_PATH at one that has it.',
            )
        endif
    elif not _rgpot_proj.get_variable('has_xtb')
        error(
            'rgpot subproject built without XTBPot. '
            + 'Reconfigure with -Drgpot:with_xtb=true, '
            + 'or wipe the build directory.',
        )
    endif
    _args += ['-DWITH_XTB']
endif

# The LAMMPS library snapshots _args. EONMPI adds a function pointer to
# LammpsLoader, and the tests compile with the final _args, so the define
# has to be set before that library() call.
if get_option('with_mpi').enabled()
    mpi_dep = dependency('mpi')
    # The server rank runs EON_SERVER_PATH through Py_Main (ClientEON.cpp),
    # so the MPI client embeds the interpreter.
    _deps += [mpi_dep, py.dependency(embed: true, required: true)]
    _args += ['-DEONMPI']
endif

# LAMMPS: always compiled, loaded at runtime via dlopen
subdir('potentials/LAMMPS')
potentials += [lammps_pot]
# features_list scans _args for -DLAMMPS_POT; the arm is unconditional.
_args += ['-DLAMMPS_POT']

if py_embed
    _args += ['-DEMBED_PYTHON']
    if get_option('with_ase')
        pyb11f_deps = [dependency('pybind11')]
        _deps += [pyb11f_deps]
        subdir('potentials/ASE')
        potentials += [ase]
        _args += ['-DWITH_ASE_POT']
    endif
endif

if get_option('with_mpi').enabled()
    # mpi and -DEONMPI were added before the LAMMPS library. This subdir
    # snapshots _deps, which already contains mpi.
    subdir('potentials/MPIPot')
    potentials += [mpipot]
endif

# Fortran-backed potentials (SW, EDIP, Lenosky, Tersoff, EAM-Al, FeHe,
# CuH2, TIP4P-H) are Fortran 2018 kernels inside librgpot, reached through
# RgpotAdapter. eOn no longer builds, installs, or dlopens Fortran itself,
# so the flang runtime handling and the per-pot .def export files retire
# with them.
if get_option('with_fortran')
    _args += ['-DWITH_FORTRAN']
endif

if get_option('with_cuh2')
    _args += ['-DCUH2_POT']
endif

if get_option('with_ase_orca')
    # TODO: Cleanup, used for Catlearn too
    # Embedding the interpreter
    pyb11f_deps = [python_dep, dependency('pybind11')]
    _deps += [pyb11f_deps]
    subdir('potentials/ASE_ORCA')
    potentials += [aseorca]
    _args += ['-DWITH_ASE_ORCA']
endif

if get_option('with_ase_nwchem')
    # TODO: Cleanup, used for Catlearn too
    # Embedding the interpreter
    python_dep = py.dependency(embed: true, required: true)
    pyb11f_deps = [
        python_dep,
        dependency('pybind11'),
        declare_dependency(link_args: '-lstdc++'),
    ]
    _deps += [pyb11f_deps]
    subdir('potentials/ASE_NWCHEM')
    potentials += [asenwchem]
    _args += ['-DWITH_ASE_NWCHEM']
endif

if get_option('with_metatomic')
    if get_option('torch_path') != ''
        # Setting up Torch..
        LIB_TORCH_PATH = get_option('torch_path')
        message('Ensure your versions are correct for metatensor and metatomic')
    elif get_option('pip_metatomic')
        # Resolve the installed torch package (site-packages layout varies by
        # venv / manylinux / no-isolation builds). Prefer import path over
        # py.get_install_dir()/torch which can miss platlib installs.
        _torch_loc = run_command(
            py,
            '-c',
            'import pathlib, torch; print(pathlib.Path(torch.__file__).resolve().parent)',
            check: true,
        )
        LIB_TORCH_PATH = _torch_loc.stdout().strip()
        message('pip_metatomic: torch at', LIB_TORCH_PATH)
    else
        LIB_TORCH_PATH = py.get_install_dir() / 'torch'
    endif
    LIB_TORCH_LIB_PATH = LIB_TORCH_PATH / 'lib'
    LIB_TORCH_INCL_PATH = LIB_TORCH_PATH / 'include'
    LIB_TORCH_INCL_PATH_API = LIB_TORCH_INCL_PATH / 'torch' / 'csrc' / 'api' / 'include'

    LIB_TORCH_INCLUDE_LIST = [
        LIB_TORCH_INCL_PATH,
        LIB_TORCH_INCL_PATH_API,
    ]


    lib_torch_list = ['c10', 'torch', 'torch_cpu']
    if not is_windows
        lib_torch_list += 'torch_global_deps'
    endif
    # On Windows, torch DLLs live in <prefix>/bin not <prefix>/lib
    LIB_TORCH_BIN_PATH = LIB_TORCH_PATH / 'bin'
    lib_torch_dirs = is_windows ? [LIB_TORCH_LIB_PATH, LIB_TORCH_BIN_PATH] : [LIB_TORCH_LIB_PATH]
    tdeps = []
    foreach lib_name : lib_torch_list
        tdeps += cppc.find_library(lib_name, dirs: lib_torch_dirs)
    endforeach

    torch_inc_args = []
    foreach d : LIB_TORCH_INCLUDE_LIST
        # Reduce the priority of each
        torch_inc_args += ['-isystem', d]
    endforeach

    torch_deps = declare_dependency(
        dependencies: tdeps,
        # Use compile arguments for -isystem
        compile_args: torch_inc_args,
    )

    if get_option('pip_metatomic')
        # Resolve site-packages + torch major.minor from the *build* interpreter
        # (py.get_install_dir() misses platlib in manylinux/no-isolation; default
        # torch_version option often lags the installed torch).
        _pip_mts = run_command(
            py,
            '-c',
            '''
import pathlib, torch
site = pathlib.Path(torch.__file__).resolve().parent.parent
tv = ".".join(torch.__version__.split(".")[:2])
# torch_version meson option override when non-default probe fails later
print(site)
print(tv)
for name in ("metatensor", "metatensor_torch", "metatomic", "vesin"):
    p = site / name
    print(name, "ok" if p.is_dir() else "missing", p)
''',
            check: true,
        )
        _pip_mts_lines = _pip_mts.stdout().strip().split('\n')
        _site_pkg = _pip_mts_lines[0].strip()
        torch_ver = _pip_mts_lines[1].strip()
        message('pip_metatomic: site-packages at', _site_pkg, 'torch_ver', torch_ver)
        # metatensor-torch >=0.10: top-level metatensor_torch/; older: metatensor/torch/
        _mts_torch_lib_new = _site_pkg / 'metatensor_torch' / f'torch-@torch_ver@' / 'lib'
        _mts_torch_lib_old = _site_pkg / 'metatensor' / 'torch' / f'torch-@torch_ver@' / 'lib'
        _mts_torch_inc_new = _site_pkg / 'metatensor_torch' / f'torch-@torch_ver@' / 'include'
        _mts_torch_inc_old = _site_pkg / 'metatensor' / 'torch' / f'torch-@torch_ver@' / 'include'
        _mts_core_lib = _site_pkg / 'metatensor' / 'lib'
        _mts_core_inc = _site_pkg / 'metatensor' / 'include'
        _mta_lib = _site_pkg / 'metatomic' / 'torch' / f'torch-@torch_ver@' / 'lib'
        _mta_inc = _site_pkg / 'metatomic' / 'torch' / f'torch-@torch_ver@' / 'include'
        _vesin_lib = _site_pkg / 'vesin' / 'lib'
        _vesin_inc = _site_pkg / 'vesin' / 'include'
        _mts_torch_lib_dirs = [_mts_torch_lib_new, _mts_torch_lib_old]
        pip_metatomic_runtime_dirs = [
            LIB_TORCH_LIB_PATH,
            _mts_core_lib,
            _mts_torch_lib_new,
            _mts_torch_lib_old,
            _mta_lib,
            _vesin_lib,
        ]
        if not is_windows
            _runtime_rpath += pip_metatomic_runtime_dirs
        endif
        mts_inc_args = []
        foreach d : [_mts_core_inc, _mts_torch_inc_new, _mts_torch_inc_old, _mta_inc]
            # MSVC has no -isystem; /external:I is opt-in. Use -I style via include_directories.
            mts_inc_args += ['-I' + d]
        endforeach
        # Windows wheels: import libs may be foo.lib, libfoo.lib, or foo.dll.lib
        # (metatensor-core win ships metatensor.dll.lib). Meson find_library('x')
        # probes x.lib / libx.lib; pass 'x.dll' to hit x.dll.lib.
        _mts_lib_candidates = is_windows ? ['metatensor', 'libmetatensor', 'metatensor.dll'] : ['metatensor']
        _mts_torch_candidates = is_windows ? ['metatensor_torch', 'libmetatensor_torch', 'metatensor_torch.dll'] : ['metatensor_torch']
        _mta_torch_candidates = is_windows ? ['metatomic_torch', 'libmetatomic_torch', 'metatomic_torch.dll'] : ['metatomic_torch']
        _mts_core_found = disabler()
        foreach _n : _mts_lib_candidates
            if not _mts_core_found.found()
                _mts_core_found = cppc.find_library(_n, dirs: [_mts_core_lib], required: false)
            endif
        endforeach
        if not _mts_core_found.found()
            error(
                'pip_metatomic: C library metatensor not found under '
                + _mts_core_lib
                + ' (tried '
                + ', '.join(_mts_lib_candidates)
                + ')',
            )
        endif
        _mts_torch_found = disabler()
        foreach _n : _mts_torch_candidates
            if not _mts_torch_found.found()
                _mts_torch_found = cppc.find_library(
                    _n,
                    dirs: _mts_torch_lib_dirs,
                    required: false,
                )
            endif
        endforeach
        if not _mts_torch_found.found()
            error(
                'pip_metatomic: C library metatensor_torch not found under '
                + ' '.join(_mts_torch_lib_dirs),
            )
        endif
        _mta_torch_found = disabler()
        foreach _n : _mta_torch_candidates
            if not _mta_torch_found.found()
                _mta_torch_found = cppc.find_library(_n, dirs: [_mta_lib], required: false)
            endif
        endforeach
        if not _mta_torch_found.found()
            error(
                'pip_metatomic: C library metatomic_torch not found under ' + _mta_lib,
            )
        endif
        mts_dep = declare_dependency(
            dependencies: [_mts_core_found, _mts_torch_found, _mta_torch_found],
            compile_args: mts_inc_args,
        )
        # vesin: single copy already on _deps via eoncbase (do not re-find pip vesin)
    else
        message(
            'Looking for existing installations of metatensor, metatensor-torch, metatomic-torch',
        )
        message(
            'If this is flaky, pass pip_metatomic=True with the torch version and let pip handle everything',
        )
        # Vesin already resolved for eoncbase (single NL backend for all pots).
        # Only metatensor / metatomic-torch are located here.
        mts_dep = dependency('metatensor_torch', required: false)
        if not mts_dep.found()
            mts_dep = dependency('metatensor-torch', required: false)
        endif
        if not mts_dep.found()
            mts_dep = dependency('metatensor_torch', method: 'cmake', required: false)
        endif
        if not mts_dep.found()
            # Manual find_library fallback
            _ms_lib = is_windows ? 'lib' : ''
            mts_dep = [
                cppc.find_library(_ms_lib + 'metatensor'),
                cppc.find_library(_ms_lib + 'metatensor_torch'),
                cppc.find_library(_ms_lib + 'metatomic_torch'),
            ]
        else
            mta_dep = dependency('metatomic_torch', required: false)
            if not mta_dep.found()
                mta_dep = dependency('metatomic-torch', required: false)
            endif
            if not mta_dep.found()
                _ms_lib = is_windows ? 'lib' : ''
                mta_dep = cppc.find_library(_ms_lib + 'metatomic_torch')
            endif
            mts_dep = [
                mts_dep,
                mta_dep,
            ]
        endif
    endif
    # Fat Metatomic (potential = Metatomic).
    # vesin already on _deps (eoncbase VesinNeighbors + classical pots).
    # mta_plugin_deps, for the RGPOT engine .so, is taken further down once
    # _deps stops growing.
    _deps += [mts_dep, torch_deps]
    _args += ['-DWITH_METATOMIC']
    metatomic_pot = library(
        'metatomic_pot',
        [
            'potentials/Metatomic/MetatomicPotential.cpp',
            'potentials/Metatomic/MetatomicCAbi.cpp',
        ],
        dependencies: _deps,
        cpp_args: _args,
        link_with: _linkto,
        include_directories: _incdirs,
        install: true,
    )
    potentials += [metatomic_pot]
endif

# NEB image forces. GCC/Clang: a bounded std::thread pool, no dependency.
# nvc++: -Dstdpar=cpu|gpu adds -stdpar=multicore|gpu and std::execution::par.
# The potential is the cost. Do not wrap Morse / LJ atom loops in
# std::for_each(std::execution::par); that is cargo-cult vs the surf-notes
# nvc++ / -mp=gpu / QMC_GPU=openmp;cuda path (QMCPACK, not eOn).
stdpar = get_option('stdpar')
if stdpar != 'none'
    cxx = meson.get_compiler('cpp')
    cxx_id = cxx.get_id()
    if cxx_id != 'nvidia_hpc' and cxx_id != 'pgi'
        error('stdpar=@0@ requires nvc++ (compiler id nvidia_hpc or pgi); got @1@'.format(stdpar, cxx_id))
    endif
    if stdpar == 'gpu'
        _args += ['-stdpar=gpu']
        _linkargs += ['-stdpar=gpu']
        # surf-notes: name the arch. QMCPACK uses QMC_GPU_ARCHS=sm_80;sm_90
        # on Snellius/Elja. nvc++ wants -gpu=cc80,cc90, not that CMake string.
        gpu_cc = get_option('stdpar_gpu_cc').strip()
        if gpu_cc != '' and gpu_cc != 'native'
            _args += ['-gpu=' + gpu_cc]
            _linkargs += ['-gpu=' + gpu_cc]
        endif
    else
        _args += ['-stdpar=multicore']
        _linkargs += ['-stdpar=multicore']
    endif
    _args += ['-DEON_PARALLEL_NEB']
endif

if get_option('with_parallel_neb') and stdpar == 'none'
    warning('with_parallel_neb no longer needs TBB and changes nothing: ' +
            'NEB images run on a bounded std::thread pool whenever ' +
            '[Main] parallel = true. Use -Dstdpar=cpu|gpu for nvc++.')
endif

if get_option('with_artn')
    # ARTn (pARTn): Fortran library for saddle point search.
    # Try system/pkg-config first, then manual artn_libdir, then cmake-build
    # the subproject. Meson's Fortran scanner breaks with -cpp + submodules,
    # so we build artn-plugin via cmake at configure time as a workaround.
    # We redeclare pARTn's C signatures inline in ARTnResource.h, so we
    # never need an -I for artn.h; hence no artn_includedir option.
    artn_dep = dependency('artn', required: false)
    if not artn_dep.found()
        artn_libdir = get_option('artn_libdir')
        artn_lib = disabler()
        if artn_libdir == ''
            # Auto-build artn-plugin from subprojects/ via cmake
            artn_src = meson.project_source_root() / 'subprojects' / 'artn-plugin'
            artn_bld = meson.project_build_root() / 'subprojects' / '_artn_cmake_build'
            message('ARTn: building subproject via cmake (meson Fortran scanner workaround)')
            run_command('cmake', '-S', artn_src, '-B', artn_bld,
                '-DCMAKE_BUILD_TYPE=Release',
                '-DWITH_LAMMPS=OFF', '-DWITH_QE=OFF', '-DWITH_SIESTA=OFF',
                check: true)
            run_command('cmake', '--build', artn_bld, '--parallel', '4', check: true)
            # LIBRARY_OUTPUT_DIRECTORY in artn-plugin's CMakeLists.txt puts
            # libartn.so at ${CMAKE_BINARY_DIR}; search a couple of plausible
            # locations symmetrically with the IRA block below in case a
            # future cmake refactor moves it.
            artn_candidates = [artn_bld, artn_bld / 'src', artn_bld / 'lib']
            foreach candidate : artn_candidates
                artn_lib = meson.get_compiler('cpp').find_library(
                    'artn',
                    dirs: [candidate],
                    required: false,
                )
                if artn_lib.found()
                    artn_libdir = candidate
                    break
                endif
            endforeach
            if not artn_lib.found()
                error('ARTn enabled but libartn was not found in the built subproject directory')
            endif
        else
            artn_lib = meson.get_compiler('cpp').find_library('artn',
                dirs: [artn_libdir], required: true)
        endif
        _runtime_rpath += [artn_libdir]
        artn_dep = declare_dependency(dependencies: artn_lib)
    endif
    _deps += [artn_dep]
    _args += ['-DWITH_ARTN']
endif

if get_option('with_ira')
    # IRA (Iterative Rotations and Assignments): Fortran library for structure comparison.
    # Same cmake workaround as ARTn.
    ira_dep = dependency('ira', required: false)
    if not ira_dep.found()
        ira_libdir = get_option('ira_libdir')
        ira_incdir = get_option('ira_includedir')
        if ira_libdir == '' or ira_incdir == ''
            ira_src = meson.project_source_root() / 'subprojects' / 'ira'
            ira_bld = meson.project_build_root() / 'subprojects' / '_ira_cmake_build'
            message('IRA: building subproject via cmake (meson Fortran scanner workaround)')
            run_command('cmake', '-S', ira_src, '-B', ira_bld,
                '-DCMAKE_BUILD_TYPE=Release',
                check: true)
            run_command('cmake', '--build', ira_bld, '--parallel', '4', check: true)
            ira_incdir = 'interface'
            ira_includes = ['-I' + ira_src / ira_incdir]
            ira_candidates = [ira_bld / 'src', ira_src / 'lib']
            ira_lib = disabler()
            foreach candidate : ira_candidates
                ira_lib = meson.get_compiler('cpp').find_library(
                    'ira',
                    dirs: [candidate],
                    required: false,
                )
                if ira_lib.found()
                    ira_libdir = candidate
                    break
                endif
            endforeach
            if not ira_lib.found()
                error('IRA enabled but libira was not found in the built or vendored library directories')
            endif
        else
            ira_includes = ['-I' + ira_incdir]
            ira_lib = meson.get_compiler('cpp').find_library(
                'ira',
                dirs: [ira_libdir],
                required: false,
            )
            if not ira_lib.found()
                error('IRA enabled but libira was not found in ira_libdir=' + ira_libdir)
            endif
        endif
        _runtime_rpath += [ira_libdir]
        ira_dep = declare_dependency(
            dependencies: ira_lib,
            compile_args: ira_includes,
        )
    endif
    _deps += [ira_dep]
    _args += ['-DWITH_IRA']
endif

if get_option('with_serve')
    # eOn as rgpot-compatible RPC *server* (clients connect to eOn; not potserv client).
    if ptlrpc_dep.found() == false
        error('with_serve requires Cap\'n Proto / rgpot ptlrpc_dep')
    endif
    _eon_build_revision = ''
    _git_exe = find_program('git', required: false)
    if _git_exe.found()
        _eon_rev = run_command(
            _git_exe,
            '-C',
            meson.project_source_root(),
            'rev-parse',
            '--short=12',
            'HEAD',
            check: false,
        )
        if _eon_rev.returncode() == 0
            _eon_build_revision = _eon_rev.stdout().strip()
        endif
    endif
    _args += [
        '-DWITH_SERVE_MODE',
        '-DEON_BUILD_VERSION="' + meson.project_version() + '"',
        '-DEON_BUILD_REVISION="' + _eon_build_revision + '"',
    ]
    eonclib_sources += ['ServeMode.cpp', 'ServeRpcServer.cpp']
    _deps += [serve_capnp_dep, ptlrpc_dep]
endif

_linkto += potentials

if get_option('with_metatomic')
    # The RGPOT metatomic engine .so compiles the same translation units as
    # eonclib, so it needs every dependency eonclib has. Snapshot here, after
    # tbb/artn/ira/capnp. readcon is already on _deps (before plugin subdirs).
    # mts_dep and torch_deps are already on _deps.
    mta_plugin_deps = _deps
endif

# version.h FEATURES must see every -DWITH_* appended by the blocks above.
_args_str = ' '.join(_args)
features_string = ''
foreach feature : features_list
    feature_flag = feature[0]
    feature_name = feature[1]
    if '-D' + feature_flag in _args_str
        features_string += '  ' + feature_name + ': enabled\n'
    else
        features_string += '  ' + feature_name + ': disabled\n'
    endif
endforeach

configure_file(
    input: 'version.h.in',
    output: 'version.h',
    configuration: {
        'VERSION': version,
        'VERSION_MAJOR': version_major,
        'VERSION_MINOR': version_minor,
        'VERSION_PATCH': version_patch,
        'BUILD_DATE': build_date,
        'OS_INFO': host_system,
        'ARCH': architecture,
        'GIT_HASH': git_hash,
        'GIT_HASH_FULL': git_hash_full,
        'FEATURES': features_string,
    },
)

# --------------------- Library

# install_rpath is the whole installed run path. Meson does not add libdir
# beside it. The client is in bindir and libeonclib is in libdir. Potential
# modules are installed into that same directory, so a library finds them
# at $ORIGIN.
_libdir = get_option('libdir')
_bindir = get_option('bindir')
if not fs.is_absolute(_libdir)
    _libdir = get_option('prefix') / _libdir
endif
if not fs.is_absolute(_bindir)
    _bindir = get_option('prefix') / _bindir
endif
_bin_to_lib = '$ORIGIN/' + fs.relative_to(_libdir, _bindir)
_lib_install_rpath = ':'.join(['$ORIGIN'] + _runtime_rpath)
_exe_install_rpath = ':'.join([_bin_to_lib] + _runtime_rpath)

eclib = library(
    'eonclib',
    sources: eonclib_sources + readcon_order_srcs,
    include_directories: _incdirs,
    dependencies: _deps + _gprd_extra_deps,
    link_with: _linkto,
    cpp_args: _args,
    link_args: _linkargs,
    build_rpath: ':'.join(_runtime_rpath),
    install_rpath: _lib_install_rpath,
    install: true,
)

# MSVC scans static libraries once, left to right. eonclib references
# the potential archives, so it has to precede them.
_linkto = [eclib] + _linkto

capnp_prog = find_program('capnp', required: true)
capnp_dep = dependency('capnp', required: true)
# pixi global exposes ~/.pixi/bin/capnp; parent(parent(that))/include
# is ~/.pixi/include and does not exist. pkg-config includedir is the
# env that actually owns /capnp/c++.capnp.
_capnp_import_path = capnp_dep.get_variable('includedir', default_value: '')
if _capnp_import_path == ''
    _capnp_import_path = fs.parent(fs.parent(capnp_prog.full_path())) / 'include'
endif
# capnp writes <schema>.c++. MSVC cl ignores that extension (D9024,
# D9027), and the librarian then fails with LNK1181 on the missing
# object. The script republishes the translation unit as .cpp.
eon_relax_capnp = custom_target(
    'eon_relax_capnp',
    input: meson.project_source_root() / 'schema' / 'eon_relax_engine.capnp',
    output: ['eon_relax_engine.capnp.h', 'eon_relax_engine.capnp.cpp'],
    command: [
        py.full_path(),
        meson.project_source_root() / 'scripts' / 'capnp_cpp_unit.py',
        '--capnp', capnp_prog.full_path(),
        '--import-path', _capnp_import_path,
        '--src-prefix', meson.project_source_root() / 'schema',
        '--schema', '@INPUT@',
        '--header', '@OUTPUT0@',
        '--source', '@OUTPUT1@',
    ],
)

# Tier-2 plugin: dlopen'ed by a learning host (gpr_optim) to run an NEB or
# a saddle search against the host's surface.
# Always a shared library. The learning host and the ABI test load it
# with dlopen/LoadLibrary. A static archive is not a loadable module:
# Windows CI sets EON_RELAX_ENGINE to libeon_relax_engine.a and
# LoadLibrary returns null. default_library=static must not apply here.
# Windows puts only readcon headers on _deps. The import library stays
# off eonclib because llvm-ar and link reject the cargo DLL. This DLL
# pulls eonclib in, so it needs the same import library the executables
# get, without the cargo sources that contain the DLL.
_engine_deps = _deps + _gprd_extra_deps + [capnp_dep]
if is_windows
    _engine_deps += readcon_dep.partial_dependency(
        includes: true,
        compile_args: true,
        link_args: true,
        links: true,
    )
endif
eon_relax_engine = shared_library(
    'eon_relax_engine',
    ['relax/RelaxEngineAbi.cpp', eon_relax_capnp],
    dependencies: _engine_deps,
    cpp_args: _args + ['-DEON_RELAX_ENGINE_BUILD'],
    link_with: _linkto,
    include_directories: _incdirs + [include_directories('.')],
    build_rpath: ':'.join(_runtime_rpath),
    install_rpath: _lib_install_rpath,
    install: true,
    gnu_symbol_visibility: 'inlineshidden',
)

# RGPOT-metatomic engine (slow path); fat pot already in potentials above.
if get_option('with_metatomic')
    metatomic_engine = library(
        'rgpot_metatomic_engine',
        [
            'potentials/Metatomic/MetatomicPotential.cpp',
            'potentials/Metatomic/MetatomicEngineAbi.cpp',
        ]
        + readcon_order_srcs,
        dependencies: mta_plugin_deps,
        cpp_args: _args + ['-DRGPOT_MTA_ENGINE_BUILD'],
        link_with: [eclib, eoncbase],
        include_directories: _incdirs,
        build_rpath: ':'.join(_runtime_rpath),
        install_rpath: _lib_install_rpath,
        install: true,
        gnu_symbol_visibility: 'default',
    )

    # One-release compat: callers and older wheels still open the previous soname.
    _engine_real = fs.name(metatomic_engine.full_path())
    _engine_compat = ['libmetatomic_engine.so']
    if host_system == 'darwin'
        _engine_compat = ['libmetatomic_engine.dylib']
    elif host_system == 'windows'
        _engine_compat = ['metatomic_engine.dll', 'libmetatomic_engine.dll']
    endif
    if host_system != 'windows'
        foreach compat : _engine_compat
            install_symlink(
                compat,
                install_dir: get_option('libdir'),
                pointing_to: _engine_real,
            )
            custom_target(
                'metatomic_engine_compat_' + compat,
                output: compat,
                command: ['ln', '-sfn', _engine_real, '@OUTPUT@'],
                depends: metatomic_engine,
                build_by_default: true,
            )
        endforeach
    endif
endif

# Full deps including readcon link_args (for executables and tests only)
_exe_deps = _deps + _gprd_extra_deps
if is_windows
    # _deps carries only readcon's partial dependency there (see the
    # partial_dependency block above); executables need the link_args too.
    # Elsewhere _deps already holds the full readcon_dep.
    #
    # Everything except sources. readcon-core declares
    # `sources: readcon_cargo`, and on Windows that array holds the DLL
    # alongside the import library. Meson hands a dependency's sources to
    # the linker, and link rejects a DLL as input:
    #   readcon_core.dll : fatal error LNK1107: invalid or corrupt file
    # link_args already names the import library by full path, which is the
    # half MSVC wants.
    _exe_deps += [
        readcon_dep.partial_dependency(
            includes: true,
            compile_args: true,
            link_args: true,
            links: true,
        ),
    ]
endif

# Export variable
eclib_dep = declare_dependency(
    link_with: [eclib, _linkto],
    include_directories: _incdirs,
    compile_args: _fargs + _args,
    dependencies: _exe_deps,
)


# --------------------- pyeonclient (nanobind)
if get_option('with_pyeonclient')
    subdir('python')
endif

# Shared eonclib already contains the static base objects. Linking that
# archive into a test as well gives macOS a second PotRegistry, so the
# process-wide force count stays at zero.
if get_option('default_library') == 'static'
    test_link = _linkto
else
    # eclib holds the one PotRegistry. Potential archives still have to be
    # on the link line: a test can construct LAMMPSPot even when the client
    # library did not pull that object in.
    test_link = [eclib] + potentials
endif

# ---------------------- Executable

ec = executable(
    'eonclient',
    sources: eonclient_sources,
    include_directories: _incdirs,
    dependencies: _exe_deps,
    link_with: _linkto,
    cpp_args: [_args, '-w'],
    link_args: _linkargs,
    build_rpath: ':'.join(_runtime_rpath),
    install_rpath: _exe_install_rpath,
    pie: false,
    # Standalone pyeonclient wheels only need the extension + shared libs.
    install: get_option('install_eon_server'),
)

if get_option('with_tests')
    if py_embed
        # ec1 = executable('testPot',
        #                  'testpot.cpp',
        #                  include_directories: _incdirs,
        #                  dependencies : _deps,
        #                  link_with : _linkto,
        #                 )
    endif

    ec2 = executable(
        'testWriteR',
        'testWriteR.cpp',
        include_directories: _incdirs,
        dependencies: _exe_deps,
        link_with: _linkto,
    )

    # An MPI client reads its flags unless an eOn server launched it, so
    # --version exits 0 outside a server and when standalone.
    if get_option('with_mpi').enabled()
        test('eonclient_mpi_version', ec, args: ['--version'], suite: ['eon'])
        test(
            'eonclient_mpi_version_standalone',
            ec,
            args: ['--version'],
            env: ['EON_CLIENT_STANDALONE=1'],
            suite: ['eon'],
        )
    endif
endif
# ------------------------ Tests

if get_option('with_tests')
    test_args = _args
    test_deps = _exe_deps
    testMain = library(
        'TestMain',
        ['unit_tests/TestMain.cpp', 'thirdparty/catch2/catch_amalgamated.cpp'],
        dependencies: _deps,  # static lib: use _deps (no readcon link_args)
        include_directories: _incdirs,
        cpp_args: test_args,
        link_with: test_link,
    )
    test_array = [  #
        ['test_impldim', 'ImpDimerTest.cpp', 'saddle_search'],
        ['test_path_integral', 'PathIntegralTest.cpp', '', 300],
        ['test_piqtst', 'PIQTSTTest.cpp', '', 120],
        ['strparse_run', 'StringHelpersTest.cpp', ''],
        ['test_matter', 'MatterTest.cpp', 'neb_morse'],
        ['test_readcon_db_mirror', 'ReadconDbMirrorTest.cpp', 'neb_morse'],
        # The cubic rate at 600 beads, and its half-ring comparison, run
        # past the default 30s.
        ['test_tunneling', 'TunnelingTest.cpp', '', 1200],
        ['test_quantum_free_energy', 'QuantumFreeEnergyTest.cpp', ''],
        ['test_pot', 'PotTest.cpp', 'neb_morse'],
        ['test_rgpot_adapter', 'RgpotAdapterTest.cpp', 'neb_morse'],
        # ['approve_params', 'ApproveParams.cpp', ''], # needs to handle additional parameters
        ['test_zbl', 'ZBLPotTest.cpp', 'au_si'],
        ['test_extpot', 'ExtPotTest.cpp', 'extpot'],
        ['test_proj_rot_trans', 'ProjectOutRotTransTest.cpp', ''],
        ['test_epicenters', 'EpiCentersTest.cpp', 'Pt_Heptamer_FrozenLayers'],
        ['test_storage_order', 'StorageOrderTest.cpp', ''],
        ['test_for_each_image', 'ForEachImageTest.cpp', ''],
        # neb_morse workdir: displacement synth test needs reactant.con (#189)
        ['test_helpers', 'HelperFunctionsTest.cpp', 'neb_morse'],
        ['test_bundling', 'BundlingTest.cpp', ''],
        ['test_xtsci_eindir', 'XtsciEindirTest.cpp', ''],
        ['test_lammps_loader', 'LammpsLoaderTest.cpp', ''],
        ['test_coverage_paths', 'CoveragePathsTest.cpp', 'neb_morse', 180],
        ['test_socket_nwchem_listen', 'SocketNWChemListenTest.cpp', ''],
        ['test_plugin_loader', 'PluginLoaderTest.cpp', ''],
        ['test_ira_resource', 'IRAResourceTest.cpp', ''],
        ['test_artn_resource', 'ARTnResourceTest.cpp', ''],
        ['test_pot_registry', 'PotRegistryTest.cpp', ''],
        ['test_runtime', 'RuntimeTest.cpp', ''],
        ['test_metatomic_loader', 'MetatomicLoaderTest.cpp', ''],
        ['test_neb_lj', 'NEBMorseTest.cpp', 'neb_morse'],
        ['test_relax_engine', 'RelaxEngineAbiTest.cpp', 'neb_morse'],
        ['test_solid_state_neb', 'SolidStateNEBTest.cpp', ''],
        ['test_neb_regression', 'NEBRegressionTest.cpp', 'neb_lj13'],
        ['test_dimer', 'DimerTest.cpp', 'neb_morse'],
        ['test_saddle', 'SaddleSearchTest.cpp', 'neb_morse'],
        ['test_confileio', 'ConFileIOTest.cpp', 'neb_morse'],
        ['test_approve_confileio', 'ApproveConFileIO.cpp', 'neb_morse'],
        ['test_minimization', 'MinimizationTest.cpp', 'neb_morse'],
        ['test_optimizer', 'OptimizerTest.cpp', ''],
        ['test_fpe_handler', 'FPEHandlerTest.cpp', ''],
        ['test_dynamics', 'DynamicsTest.cpp', 'neb_morse'],
        # The LJ13 rate is three cooled ring searches.
        ['test_jobs', 'JobIntegrationTest.cpp', 'neb_morse', 1200],
        ['test_ohtst_scan', 'OHTSTScanTest.cpp', ''],
        ['test_hessian', 'HessianTest.cpp', 'neb_morse'],
        ['test_geom', 'GeometryAnalysisTest.cpp', 'neb_morse'],
        ['test_params_json', 'ParametersJSONTest.cpp', ''],
        ['test_params_ssot', 'ParametersSSOTTest.cpp', ''],
        ['test_bondboost', 'BondBoostTest.cpp', 'neb_morse'],
        ['test_job_result', 'JobResultProvenanceTest.cpp', ''],
        ['test_prefactor', 'PrefactorTest.cpp', 'neb_morse'],
        ['test_job_runners', 'JobRunnerCoverageTest.cpp', 'neb_morse'],
        ['test_emt_cu', 'EMTCuTest.cpp', 'cu_fcc'],
        ['test_neb_subsystem', 'NEBSubsystemTest.cpp', ''],
        ['test_neb_zoom', 'NEBZoomTest.cpp', ''],
    ]
    # posix_spawn. On Windows the only case is skipped and Catch2 exits 4.
    if not is_windows
        test_array += [['test_client_timing', 'ClientTimingTest.cpp', '', 60]]
    endif
    # Fortran potential tests: require Fortran shared libs at runtime
    if get_option('with_fortran')
        test_array += [
            ['test_si_pot', 'SiPotTest.cpp', 'si_diamond'],
            ['test_eam_al', 'EAMAlTest.cpp', 'al_fcc'],
            ['test_fehe', 'FeHeTest.cpp', 'fe_bcc'],
        ]
    endif
    if get_option('with_ase')
        test_array += [['test_ase_pot', 'ASEPotTest.cpp', 'ase_pot']]
    endif
    if get_option('with_metatomic')
        test_array += [['test_mta', 'MetatomicTest.cpp', 'lj38']]
        # DynLib.h (LoadLibrary / dlopen) — builds on MSVC and POSIX.
        test_array += [['test_mta_engine_abi', 'MetatomicEngineAbiTest.cpp', '']]
    endif
    if get_option('with_xtb')
        test_array += [['test_xtb', 'XTBTest.cpp', 'sulfolene']]
        test_array += [['test_cineb_xtb', 'CINEBXTBTest.cpp', 'cineb_xtb', 120]]
    endif
    if get_option('with_serve')
        test_array += [['test_serve_spec', 'ServeSpecParseTest.cpp', '']]
        test_array += [['test_rpc_capabilities', 'RpcCapabilitiesTest.cpp', '']]
    endif
    if get_option('with_artn')
        test_array += [['test_artn', 'ARTnTest.cpp', 'neb_morse']]
    endif
    if get_option('with_rgsaddle')
        test_array += [['test_rgsaddle_band', 'RgsaddleBandTest.cpp', 'neb_morse']]
    endif
    if get_option('with_ira')
        test_array += [['test_ira', 'IRATest.cpp', 'saddle_search']]
    endif
    if with_gprd
        test_array += [['test_gpr_helpers', 'GPRHelpersTest.cpp', '']]
        # Frozen GprMeanPotential fixed-point parity between the eonclient
        # LBFGS and the gpr band engine.
        test_array += [['test_gprmean_fixedpoint', 'GprMeanFixedPointTest.cpp', '']]
    endif
    # Build the Fortran potential library search path for tests.
    # Each Fortran shared lib lives in its subdir of the build tree,
    # matching the source layout: build/client/potentials/<Name>/
    _fortran_pot_subdirs = ['Aluminum', 'EDIP', 'FeHe', 'Lenosky', 'SW', 'Tersoff']
    _fortran_pot_dirs = []
    if get_option('with_fortran')
        foreach sub : _fortran_pot_subdirs
            _fortran_pot_dirs += meson.project_build_root() / 'client' / 'potentials' / sub
        endforeach
    endif
    # Metatomic C-ABI plugin lives next to MetatomicPotential sources in the build tree.
    if get_option('with_metatomic')
        _fortran_pot_dirs += meson.project_build_root() / 'client'
        _fortran_pot_dirs += meson.project_build_root() / 'client' / 'potentials' / 'Metatomic'
    endif
    _test_env = environment()
    _systems_dir = join_paths(
        meson.project_source_root(),
        'client',
        'unit_tests',
        'data',
        'systems',
    )
    _test_env.set('EON_TEST_SYSTEMS_DIR', _systems_dir)
    _test_env.set('EONCLIENT', ec.full_path())
    _test_env.set('EON_RELAX_ENGINE', eon_relax_engine.full_path())
    # The Metatomic job puts torch/lib ahead of the pixi prefix. The
    # engine DLL was linked against the prefix's capnp, kj, and flang
    # runtimes. A same-named DLL from torch is missing those entry
    # points and the process exits 0xC0000139. The prefix bins, then
    # the built readcon DLL, come first.
    if is_windows
        _prefix_bin = get_option('prefix') / 'bin'
        _prefix_lbin = get_option('prefix') / 'Library' / 'bin'
        if fs.is_dir(_prefix_bin)
            _test_env.prepend('PATH', _prefix_bin)
        endif
        if fs.is_dir(_prefix_lbin)
            _test_env.prepend('PATH', _prefix_lbin)
        endif
    endif
    if _readcon_runtime_dir != ''
        _test_env.prepend('PATH', _readcon_runtime_dir)
    endif
    # Windows loads a DLL's imports from the executable directory before
    # PATH. Stage the engine's non-system imports beside the test, the
    # same place the Fortran pot DLLs sit. Torch on PATH cannot shadow
    # them.
    _stage_engine_imports = []
    _have_staged_engine_imports = false
    if is_windows
        _stage_search = []
        if _readcon_runtime_dir != ''
            _stage_search += _readcon_runtime_dir
        endif
        _lbin = get_option('prefix') / 'Library' / 'bin'
        _pbin = get_option('prefix') / 'bin'
        if fs.is_dir(_lbin)
            _stage_search += _lbin
        endif
        if fs.is_dir(_pbin)
            _stage_search += _pbin
        endif
        _stage_engine_imports = custom_target(
            'stage-engine-dll-imports',
            input: eon_relax_engine,
            output: 'engine-dll-imports.stamp',
            command: [
                py,
                meson.project_source_root() / 'scripts' / 'ci' / 'win_stage_imports.py',
                '@INPUT@',
                '@OUTPUT@',
            ] + _stage_search,
            build_by_default: true,
        )
        _have_staged_engine_imports = true
    endif
    if _fortran_pot_dirs.length() > 0
        # FortranPotLoader splits on ';' on Windows (drive letters break ':').
        # Same pathsep rule as PATH: win=';' unix=':'.
        _pot_pathsep = is_windows ? ';' : ':'
        _test_env.set(
            'EON_POTENTIALS_PATH',
            _pot_pathsep.join(_fortran_pot_dirs),
        )
    endif

    foreach test : test_array
        test_timeout = test.get(3, 30)
        _case_env = _test_env
        _case_depends = []
        if _have_readcon_dll_beside
            _case_depends += [_readcon_dll_beside]
        endif
        if _have_staged_engine_imports
            _case_depends += [_stage_engine_imports]
        endif
        if test.get(0) == 'test_client_timing'
            _case_depends += [ec]
        elif test.get(0) == 'test_serve_spec'
            # Spawn eonclient. The shared test env is one object, so this
            # case gets its own copy plus the client path.
            _case_env = environment()
            _case_env.set('EON_TEST_SYSTEMS_DIR', _systems_dir)
            _case_env.set('EONCLIENT', ec.full_path())
            _case_depends = [ec]
        endif
        test(
            test.get(0),
            executable(
                test.get(0),
                sources: [
                    'unit_tests/' + test.get(1),
                    'thirdparty/catch2/catch_amalgamated.cpp',
                ],
                include_directories: test.get(0) == 'test_relax_engine'
                ? _incdirs + [include_directories('.')]
                : _incdirs,
                cpp_args: test_args,
                # Include /STACK:16M (Windows) so EAM Al Fortran pots do not
                # SIGSEGV with STATUS_STACK_OVERFLOW (default stack is 1 MB).
                link_args: _linkargs,
                dependencies: test.get(0) == 'test_relax_engine'
                ? test_deps + [capnp_dep]
                : test_deps,
                link_with: test.get(0) == 'test_relax_engine'
                ? test_link + [eon_relax_engine]
                : test_link,
                build_rpath: ':'.join(_runtime_rpath),
            ),
            workdir: meson.project_source_root() + '/client/unit_tests/data/systems/' + test.get(
                2,
            ),
            timeout: test_timeout,
            env: _case_env,
            depends: _case_depends,
        )
    endforeach
    if get_option('with_water')
        test(
            'test_water_cutover',
            executable(
                'test_water_cutover',
                sources: [
                    'unit_tests/WaterCutoverTest.cpp',
                    'thirdparty/catch2/catch_amalgamated.cpp',
                ],
                include_directories: _incdirs,
                cpp_args: test_args,
                link_args: _linkargs,
                dependencies: test_deps,
                link_with: test_link,
                build_rpath: ':'.join(_runtime_rpath),
            ),
            timeout: 60,
            env: _test_env,
        )
    endif
    if get_option('with_xtb')
        # Timing harness for linked ship vs rgpot dlopen (not a Catch suite).
        executable(
            'xtb_ship_bench',
            sources: ['unit_tests/xtb_ship_bench.cpp'],
            dependencies: test_deps,
            include_directories: _incdirs,
            cpp_args: test_args,
            link_with: test_link,
            build_rpath: ':'.join(_runtime_rpath),
            install: false,
        )
    endif
    # Special handling for the nwchem test (due to the client server bit)
    nwchem_exe = find_program('nwchem', required: false)
    bash_exe = find_program('bash', required: false)
    # --- Special handling for the SocketNWChem integration test ---
    if nwchem_exe.found()
        test_socket_nwchem_exe = executable(
            'test_socket_nwchem',
            sources: [
                'unit_tests/SocketNWChemPotTest.cpp',
                'thirdparty/catch2/catch_amalgamated.cpp',
            ],
            dependencies: test_deps,
            include_directories: _incdirs,
            cpp_args: test_args,
            link_with: test_link,
        )

        test_workdir = meson.project_source_root() + '/client/unit_tests/data/systems/nwchem_test'

        test(
            'test_socket_nwchem',
            bash_exe,
            args: [
                meson.project_source_root() + '/client/unit_tests/run_nwchem_test.sh',
                ec,
                nwchem_exe,
                'nwchem_socket.nwi',  # The input file for the client
                test_workdir,
            ],
            # Note: The workdir is now an argument to the script, not a meson property for the test
            timeout: 120,
        )
    else
        message('NWChem executable not found, skipping SocketNWChemPotTest.')
    endif

    if not is_windows
        test_rgpot_exe = executable(
            'test_rgpot_pot',
            sources: [
                'unit_tests/RgpotPotTest.cpp',
                'unit_tests/CpmdMessageTest.cpp',
                'unit_tests/GroupScheduleTest.cpp',
                'thirdparty/catch2/catch_amalgamated.cpp',
            ],
            dependencies: test_deps,
            include_directories: _incdirs,
            cpp_args: test_args,
            link_with: test_link,
            build_rpath: ':'.join(_runtime_rpath),
        )
        _rgpot_test_env = environment()
        _rgpot_test_env.set(
            'EON_TEST_SYSTEMS_DIR',
            join_paths(
                meson.project_source_root(),
                'client',
                'unit_tests',
                'data',
                'systems',
            ),
        )
        # Optional nwchemc/cpmdc engines: tests SKIP when absent. Catch2 exits 4
        # when every case is skipped unless --allow-running-no-tests is set.
        test(
            'test_rgpot_pot',
            test_rgpot_exe,
            args: ['--allow-running-no-tests'],
            timeout: 120,
            suite: 'rgpot',
            is_parallel: false,
            env: _rgpot_test_env,
            workdir: meson.project_source_root() +
                     '/client/unit_tests/data/systems/rgpot_nwchem_test',
        )
        # librgpot and librgpot_pot are MPI-free; the MPI side loads on
        # demand, so eonclient itself must not map libmpi or a transport.
        # eOn's own with_mpi links eonclient to libmpi, so that check
        # holds only when that option is off.
        if (
            _rgpot_mpi_dep.found()
            and host_machine.system() == 'linux'
            and not get_option('with_mpi').enabled()
        )
            test(
                'test_eonclient_no_mpi_at_load',
                find_program('bash'),
                args: [
                    meson.project_source_root()
                    / 'client'
                    / 'unit_tests'
                    / 'check_no_mpi_closure.sh',
                    ec.full_path(),
                ],
                depends: [ec],
                suite: ['eon', 'rgpot'],
            )
        endif
        # The three cases call the calculator-group agreement. That code
        # is compiled only when rgpot links MPI. mpirun on PATH is not
        # that build, and the cases then fail closed.
        _rgpot_mpirun = find_program('mpirun', 'mpiexec', required: false)
        _rgpot_bash = find_program('bash', required: false)
        if _rgpot_mpi_dep.found() and _rgpot_mpirun.found() and _rgpot_bash.found()
            rgpot_mpi_exit = executable(
                'rgpot_mpi_exit',
                sources: ['unit_tests/rgpot_mpi_exit.cpp'],
                dependencies: test_deps,
                include_directories: _incdirs,
                cpp_args: test_args,
                link_args: _linkargs,
                link_with: test_link,
                build_rpath: ':'.join(_runtime_rpath),
            )
            _rgpot_mpi_script = (
                meson.project_source_root()
                / 'client'
                / 'unit_tests'
                / 'run_rgpot_mpi.sh'
            )
            test(
                'test_rgpot_mpi_fault',
                _rgpot_bash,
                args: [_rgpot_mpi_script, rgpot_mpi_exit.full_path(), 'fault'],
                timeout: 90,
                suite: ['eon', 'rgpot'],
                is_parallel: false,
                env: _rgpot_test_env,
            )
            test(
                'test_rgpot_mpi_abort',
                _rgpot_bash,
                args: [_rgpot_mpi_script, rgpot_mpi_exit.full_path(), 'abort'],
                timeout: 90,
                suite: ['eon', 'rgpot'],
                is_parallel: false,
                env: _rgpot_test_env,
            )
            test(
                'test_rgpot_mpi_single',
                _rgpot_bash,
                args: [_rgpot_mpi_script, rgpot_mpi_exit.full_path(), 'single'],
                timeout: 90,
                suite: ['eon', 'rgpot'],
                is_parallel: false,
                env: _rgpot_test_env,
            )
            test(
                'test_rgpot_mpi_uneven',
                _rgpot_bash,
                args: [_rgpot_mpi_script, rgpot_mpi_exit.full_path(), 'uneven'],
                timeout: 90,
                suite: ['eon', 'rgpot'],
                is_parallel: false,
                env: _rgpot_test_env,
            )
            test(
                'test_rgpot_mpi_params',
                _rgpot_bash,
                args: [_rgpot_mpi_script, rgpot_mpi_exit.full_path(), 'params'],
                timeout: 90,
                suite: ['eon', 'rgpot'],
                is_parallel: false,
                env: _rgpot_test_env,
            )
        endif
    endif

    # Benchmarks (not in default test suite -- run via: meson test --benchmark)
    bench_pot_exe = executable('bench_pot',
        [
            'unit_tests/PotBenchmark.cpp',
            'unit_tests/TestMain.cpp',
            'thirdparty/catch2/catch_amalgamated.cpp',
        ],
        dependencies: test_deps,
        include_directories: _incdirs,
        cpp_args: test_args,
        link_with: test_link,
    )
    benchmark('bench_pot', bench_pot_exe,
        args: ['[.benchmark]'],
        workdir: meson.project_source_root() + '/client/unit_tests/data/systems/neb_morse',
        timeout: 120,
    )
endif
