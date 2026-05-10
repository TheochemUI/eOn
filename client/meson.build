_incdirs += [include_directories('thirdparty')]
_runtime_rpath = []

cmake = import('cmake')

# Argum (C++20 argument parser) - single header in thirdparty

add_languages('c', required: true)
cc = meson.get_compiler('c')
cppc = meson.get_compiler('cpp')

use_fortran = false
if get_option('with_fortran') or get_option('with_cuh2')
    use_fortran = true
    add_languages('fortran', required: true)
    fc = meson.get_compiler('fortran')
    _fargs += fc.get_supported_arguments(['-fno-implicit-none'])
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

# TODO(rg): use fc.get_supported_arguments and add more compiler support
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
            '-w',
            # ^--- Silence.. TODO(rg): Seriously a very bad idea long term
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

# Release optimization flags
if get_option('buildtype') == 'release' or get_option('buildtype') == 'debugoptimized'
    _opt_args = cppc.get_supported_arguments(['-O3', '-flto=auto'])
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
version_output = run_command(version_script).stdout().strip()
version_bits = version_output.split(',')
version = version_bits[0]
build_date = version_bits[1]
git_hash = version_bits[2]

# Parse semver into components
version_parts = version.split('.')
version_major = version_parts[0].to_int()
version_minor = version_parts[1].to_int()
version_patch = version_parts[2].split('+')[0].to_int()  # Remove any +git suffix

architecture = host_machine.cpu_family()
message(host_system)
message(architecture)

# Generate features string for version.h.in.
#
# Each entry is (compile_flag, display_name, mode) where mode is one of:
#   'compile'    -- check '-Dflag' in _args (compile-time CPP define).
#                   Reports 'enabled' (compile-linked) or 'disabled'.
#   'runtime'    -- always compiled into the client; loaded via dlopen
#                   at run time (LAMMPS, ARTn, IRA, XTB). Always reports
#                   'enabled (dlopen at runtime)' since the shim is
#                   present, irrespective of whether the dynamic library
#                   is found at run time. Call sites use is_loaded() /
#                   require_loaded() to gate the runtime path.
#   'subprocess' -- POSIX-only subprocess shim (VASP). Always compiled
#                   on POSIX, omitted on Windows.
features_list = [
    ['WITH_GPRD',          'GPRD',           'compile'],
    ['WITH_GP_SURROGATE',  'GP Surrogate',   'compile'],
    ['WITH_CATLEARN',      'CatLearn',       'compile'],
    # VASP is unconditionally compiled on POSIX (omitted on Windows).
    # The class spawns ./runvasp.sh and reads FU/POSCAR files; nothing
    # is linked or dlopen'd. The banner reports the shim's presence,
    # not whether VASP itself is on PATH.
    ['VASP_POT',           'VASP',           'subprocess'],
    ['WITH_WATER',         'Water',          'compile'],
    ['WITH_AMS',           'AMS',            'compile'],
    # XTB: always compiled, dlopen at runtime via XtbLoader.
    ['XTB_POT',            'XTB',            'runtime'],
    # LAMMPS is unconditionally compiled (see subdir('potentials/LAMMPS')
    # below). The dlopen lookup of liblammps.so happens at run time inside
    # LammpsLoader::instance(); the banner just tells the user the shim is
    # present, not whether liblammps was located.
    ['LAMMPS_POT',         'LAMMPS',         'runtime'],
    ['WITH_ASE_POT',       'ASE Potentials', 'compile'],
    ['EONMPI',             'MPI',            'compile'],
    ['WITH_FORTRAN',       'Fortran',        'compile'],
    ['CUH2_POT',           'CuH2',           'compile'],
    ['WITH_ASE_ORCA',      'ASE-ORCA',       'compile'],
    ['WITH_ASE_NWCHEM',    'ASE-NWChem',     'compile'],
    # Metatomic is compile-linked (subdir('potentials/Metatomic') links
    # the metatomic_pot library directly via library() with link_with).
    # No dlopen; static check on -DWITH_METATOMIC is correct.
    ['WITH_METATOMIC',     'Metatomic',      'compile'],
    ['WITH_SERVE_MODE',    'Serve Mode',     'compile'],
    ['EMBED_PYTHON',       'Embedded Python','compile'],
    # ARTn / IRA: same pattern as LAMMPS. The *Resource shim is
    # unconditionally compiled; libartn.so / libira.so are dlopen'd
    # inside *Resource::instance(); is_loaded() / require_loaded()
    # gate the call sites at run time.
    ['ARTN_POT',           'ARTn',           'runtime'],
    ['IRA_POT',            'IRA',            'runtime'],
]

# Join _args into a single string for checking
_args_str = ' '.join(_args)

features_string = ''
foreach feature : features_list
    feature_flag = feature[0]
    feature_name = feature[1]
    feature_mode = feature.length() > 2 ? feature[2] : 'compile'
    if feature_mode == 'runtime'
        features_string += '  ' + feature_name + ': enabled (dlopen at runtime)\n'
    elif feature_mode == 'subprocess'
        # POSIX-only subprocess shims (VASP). Always compiled where the
        # platform supports it; the user supplies the external binary
        # at run time on PATH.
        if host_system == 'windows'
            features_string += '  ' + feature_name + ': disabled (POSIX only)\n'
        else
            features_string += '  ' + feature_name + ': enabled (subprocess)\n'
        endif
    elif '-D' + feature_flag in _args_str
        features_string += '  ' + feature_name + ': enabled\n'
    else
        features_string += '  ' + feature_name + ': disabled\n'
    endif
endforeach

# Configure the header file
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
        'FEATURES': features_string,
    },
)

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

quill_dep = dependency('quill', required: true)
_deps += quill_dep
_args += ['-DQUILL_DISABLE_NON_PREFIXED_MACROS']

inireader_dep = dependency('INIReader', required: false)
if not inireader_dep.found()
    inireader_dep = dependency(
        'inireader',
        fallback: ['inih', 'INIReader_dep'],
        required: true,
        default_options: ['default_library=static'],
    )
endif
_deps += inireader_dep

# Always use the subproject: system nlohmann_json cmake configs can export
# -I/usr/include which breaks conda-forge sysroot isolation on rolling distros.
nlohmann_json_proj = subproject(
    'nlohmann_json',
    default_options: ['default_library=static'],
)
nlohmann_json_dep = nlohmann_json_proj.get_variable('nlohmann_json_dep')
_deps += nlohmann_json_dep

# Highway SIMD (optional, for vectorized potentials)
hwy_dep = dependency('libhwy', required: false)
if not hwy_dep.found() and host_machine.system() != 'windows'
    # cmake subproject produces .a archives that MSVC link.exe cannot read
    message('Highway: system libhwy not found, trying cmake subproject...')
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
    endif
endif
if hwy_dep.found()
    _deps += hwy_dep
    _args += ['-DWITH_HIGHWAY']
    message('Highway SIMD: enabled')
else
    message('Highway SIMD: disabled (not found)')
endif

eoncbase_sources = ['fpe_handler.cpp', 'PotRegistry.cpp']

eoncbase = library(
    'eoncbase',
    sources: eoncbase_sources,
    dependencies: _deps,
    cpp_args: _args,
    include_directories: _incdirs,
    install: true,
)

_linkto += eoncbase

# Potentials
subdir('potentials/EAM')
subdir('potentials/EMT')
subdir('potentials/ExtPot')
subdir('potentials/IMD')
subdir('potentials/LJ')
subdir('potentials/LJCluster')
subdir('potentials/Morse')
subdir('potentials/ZBL')
if not is_windows
    subdir('potentials/SocketNWChem')
endif
subdir('potentials/NewPot')

potentials = [
    eam,
    emt,
    extpot,
    imd,
    lennard_jones,
    lennard_jones_cluster,
    morse,
    newpot,
    zbl,
]
if not is_windows
    potentials += [socket_nwchem]
endif

eonclib_sources = [
    'Parameters.cpp',
    'ParametersINI.cpp',
    'ParametersJSON.cpp',
    'ConFileIO.cpp',
    'GeometryAnalysis.cpp',
    'Optimizer.cpp',
    'IDPPObjectiveFunction.cpp',
    'PrefactorJob.cpp',
    'LBFGS.cpp',
    'ReplicaExchangeJob.cpp',
    'BondBoost.cpp',
    'Job.cpp',
    'GlobalOptimization.cpp',
    'LowestEigenmode.cpp',
    'MinModeSaddleSearch.cpp',
    'StructureComparisonJob.cpp',
    'SteepestDescent.cpp',
    'ImprovedDimer.cpp',
    'PointJob.cpp',
    'Prefactor.cpp',
    'ConjugateGradients.cpp',
    'Matter.cpp',
    'FiniteDifferenceJob.cpp',
    'Lanczos.cpp',
    'HessianJob.cpp',
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
    'NEBSpringForce.cpp',
    'NEBTangent.cpp',
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
    'MatrixHelpers.hpp',  # Template
    'Dimer.cpp',
    'Dynamics.cpp',
    'GlobalOptimizationJob.cpp',
    'BiasedGradientSquaredDescent.cpp',
    'SafeHyperJob.cpp',
    'MinimizationJob.cpp',
    'Quickmin.cpp',
    'ParallelReplicaJob.cpp',
    'Potential.cpp',
    'SurrogatePotential.cpp',  # Part of the interface
    'BasinHoppingJob.cpp',
    'FIRE.cpp',
    'EpiCenters.cpp',
    'SaddleSearchJob.cpp',
    'BasinHoppingSaddleSearch.cpp',
]

eonclient_sources = ['ClientEON.cpp', 'CommandLine.cpp']

# ------------------- Options

if get_option('with_gprd')
    gprd_opts = {
        'default_library': 'static',
        'use_mkl': get_option('use_mkl'),
        'use_hdf5': true,
        'do_install': true,
        'buildtype': 'release',
        'linalg_backend': get_option('gprd_linalg_backend'),
    }
    libgprd_proj = subproject('gpr_optim', default_options: gprd_opts)
    libgprd = libgprd_proj.get_variable('gpr_optim_dep')
    gprd_deps = libgprd_proj.get_variable('gpr_optim_deps')
    # subdir('potentials/GPRPotential')
    # potentials += [ gprpot ]
    _args += ['-DWITH_GPRD']
    eonclib_sources += ['AtomicGPDimer.cpp', 'GPRHelpers.cpp']
    _deps += [libgprd, gprd_deps]
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

# VASP: always compiled on POSIX (Windows omits via #ifndef _WIN32 in
# Potential.cpp). The shim spawns `runvasp.sh` at use time; no library
# dependency at build time. The user gets VASP support iff the script is
# on PATH and PotType::VASP is requested.
if host_system != 'windows'
    subdir('potentials/VASP')
    potentials += [vasp]
endif

if get_option('with_water')
    subdir('potentials/Water')
    subdir('potentials/Water_Pt')
    potentials += [water, water_pt]
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

# XTB: always compiled, dlopens libxtb at runtime via XtbLoader.
# Mirrors the LAMMPS / ARTn / IRA pattern; no link-time xtb dependency.
subdir('potentials/XTBPot')
potentials += xtb_eon

# LAMMPS: always compiled, loaded at runtime via dlopen
subdir('potentials/LAMMPS')
potentials += [lammps_pot]

if py_embed
    _args += ['-DEMBED_PYTHON']
    # subdir('potentials/PyAMFF')
    # pyamff,
    # subdir('potentials/QSC')
    # potentials += [ qsc ]
    if get_option('with_ase')
        pyb11f_deps = [dependency('pybind11')]
        _deps += [pyb11f_deps]
        subdir('potentials/ASE')
        potentials += [ase]
        _args += ['-DWITH_ASE_POT']
    endif
endif

if get_option('with_mpi')
    subdir('potentials/MPIPot')
    potentials += [mpipot]
    _args += ['-DEONMPI']
    mpi_dep = dependency('mpi')
    _deps += [mpi_dep]
endif

if get_option('with_fortran')
    subdir('potentials/Aluminum')
    subdir('potentials/EDIP')
    subdir('potentials/FeHe')
    subdir('potentials/Lenosky')
    subdir('potentials/SW')
    subdir('potentials/Tersoff')
    if get_option('with_water')
        subdir('potentials/Water_H')
        potentials += [water_h]
    endif
    potentials += [aluminum, edip, fehe, lenosky, sw, tersoff]
    _args += ['-DWITH_FORTRAN']
endif

if get_option('with_cuh2')
    subdir('potentials/CuH2')
    potentials += [cuh2]
    _args += ['-DCUH2_POT', '-DWITH_FORTRAN']
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
        torch_ver = get_option('torch_version')
        pip_metatomic_runtime_dirs = [
            LIB_TORCH_LIB_PATH,
            py.get_install_dir() / 'metatensor' / 'lib',
            py.get_install_dir() / 'metatensor' / 'torch' / f'torch-@torch_ver@' / 'lib',
            py.get_install_dir() / 'metatomic' / 'torch' / f'torch-@torch_ver@' / 'lib',
            py.get_install_dir() / 'vesin' / 'lib',
        ]
        if not is_windows
            _runtime_rpath += pip_metatomic_runtime_dirs
        endif
        mts_inc_args = []
        foreach d : [
            py.get_install_dir() / 'metatensor' / 'include',
            py.get_install_dir() / 'metatensor' / 'torch' / f'torch-@torch_ver@' / 'include',
            py.get_install_dir() / 'metatomic' / 'torch' / f'torch-@torch_ver@' / 'include',
        ]
            mts_inc_args += ['-isystem', d]
        endforeach
        mts_dep = declare_dependency(
            dependencies: [
                cppc.find_library(
                    'metatensor',
                    dirs: [py.get_install_dir() / 'metatensor' / 'lib'],
                ),
                cppc.find_library(
                    'metatensor_torch',
                    dirs: [
                        py.get_install_dir() / 'metatensor' / 'torch' / f'torch-@torch_ver@' / 'lib',
                    ],
                ),
                cppc.find_library(
                    'metatomic_torch',
                    dirs: [
                        py.get_install_dir() / 'metatomic' / 'torch' / f'torch-@torch_ver@' / 'lib',
                    ],
                ),
            ],
            compile_args: mts_inc_args,
        )
        vesin_dep = cppc.find_library(
            'vesin',
            dirs: [py.get_install_dir() / 'vesin' / 'lib'],
            required: false,
        )
        if vesin_dep.found()
            vesin_dep = declare_dependency(
                dependencies: [
                    cppc.find_library(
                        'vesin',
                        dirs: [py.get_install_dir() / 'vesin' / 'lib'],
                        required: false,
                    ),
                ],
                compile_args: ['-isystem', py.get_install_dir() / 'vesin' / 'include'],
            )
        endif
    else
        message(
            'Looking for existing installations of metatensor, metatensor-torch, metatomic-torch and vesin',
        )
        message(
            'If this is flaky, pass pip_metatomic=True with the torch version and let pip handle everything',
        )
        # Try system dependency() first (pkg-config / cmake), fall back to find_library
        # xref: https://github.com/metatensor/metatensor/issues/930
        # On Windows (MSVC), static libraries are often prefixed with 'lib'
        vesin_dep = dependency('vesin', required: false)
        if not vesin_dep.found()
            vesin_dep = cppc.find_library('vesin', required: false)
        endif
        if not vesin_dep.found() and is_windows
            vesin_dep = cppc.find_library('libvesin', required: false)
        endif
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
    if not vesin_dep.found()
        vesin_int = library(
            'vesin_internal',
            ['thirdparty/vesin/vesin-single-build.cpp'],
            install: true,
        )
        vesin_dep = declare_dependency(
            link_with: vesin_int,
            include_directories: include_directories('thirdparty/vesin'),
        )
    endif
    _deps += [vesin_dep, mts_dep, torch_deps]
    subdir('potentials/Metatomic')
    potentials += [metatomic_pot]
    _args += ['-DWITH_METATOMIC']
endif

if get_option('with_parallel_neb')
    tbb_dep = dependency('tbb', required : true)
    _deps += [tbb_dep]
    _args += ['-DEON_PARALLEL_NEB']
endif

# ARTn (pARTn) and IRA: always compiled, dlopen'd at runtime.
#
# The *Resource shims (libs/ARTn/ARTnResource.cpp, libs/IRA/IRAResource.cpp)
# use dlopen(libartn.so) / dlopen(libira.so) when require_loaded() is first
# called. Both shims redeclare the C entry points inline (see
# ARTnResource.h / IRAResource.h), so neither needs an external header
# at compile time. Users supply libartn.so / libira.so on
# LD_LIBRARY_PATH (or DYLD_LIBRARY_PATH on macOS) and ARTnSaddleSearch /
# IRACompare check is_loaded() at the call site -- identical pattern to
# LAMMPS today.
eonclib_sources += [
    'ARTnSaddleSearch.cpp',
    'libs/ARTn/ARTnResource.cpp',
    'IRACompare.cpp',
    'libs/IRA/IRAResource.cpp',
]

if get_option('with_serve')
    # rgpot-compatible RPC serve mode: wraps any eOn potential and serves it
    # over Cap'n Proto RPC using the rgpot Potentials.capnp schema.
    # Only links capnp-rpc and the generated schema (ptlrpc); does NOT link
    # rgpot types to avoid the AtomMatrix name collision with eOn's Eigen.h.
    serve_capnp_dep = dependency('capnp-rpc', required: true)
    rgpot_serve_opts = {
        'with_rpc_client_only': true,
        'pure_lib': true,
        'with_tests': false,
        'with_examples': false,
    }
    rgpot_proj = subproject('rgpot', default_options: rgpot_serve_opts)
    ptlrpc_dep = rgpot_proj.get_variable('ptlrpc_dep')
    _args += ['-DWITH_SERVE_MODE']
    eonclib_sources += ['ServeMode.cpp', 'ServeRpcServer.cpp']
    _deps += [serve_capnp_dep, ptlrpc_dep]
endif

_linkto += potentials

# readcon-core: Rust-based .con file I/O (required)
readcon_proj = subproject(
    'readcon-core',
    default_options: {
        'with_tests': false,
        'with_examples': false,
        'default_library': 'static',
    },
)
readcon_dep = readcon_proj.get_variable('readcon_dep')
if is_windows
    # On Windows, readcon_dep's sources (custom_target .lib) get passed to
    # llvm-ar which misinterprets the .lib as the archive output name.
    # Strip sources and link_args for the static library; the full dep
    # with link_args goes to executables via _exe_deps.
    # Build ordering is ensured by the `depends` kwarg on the library target.
    readcon_compile_dep = readcon_dep.partial_dependency(
        includes: true,
        compile_args: true,
    )
    _deps += [readcon_compile_dep]
    # No need for explicit depends: cbindgen in-tree header already exists
else
    _deps += [readcon_dep]
endif

# --------------------- Library

eclib = library(
    'eonclib',
    sources: eonclib_sources,
    include_directories: _incdirs,
    dependencies: _deps,
    link_with: _linkto,
    cpp_args: _args,
    link_args: _linkargs,
    build_rpath: ':'.join(_runtime_rpath),
    install_rpath: ':'.join(_runtime_rpath),
    install: true,
)

_linkto += eclib

# Full deps including readcon link_args (for executables and tests only)
_exe_deps = _deps + [readcon_dep]

# Export variable
eclib_dep = declare_dependency(
    link_with: [eclib, _linkto],
    include_directories: _incdirs,
    compile_args: _fargs + _args,
    dependencies: _exe_deps,
)

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
    install_rpath: ':'.join(_runtime_rpath),
    install: true,
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
        link_with: _linkto,
    )
    test_array = [  #
        ['test_impldim', 'ImpDimerTest.cpp', 'saddle_search'],
        ['strparse_run', 'StringHelpersTest.cpp', ''],
        ['test_matter', 'MatterTest.cpp', 'neb_morse'],
        ['test_pot', 'PotTest.cpp', 'neb_morse'],
        # ['approve_params', 'ApproveParams.cpp', ''], # needs to handle additional parameters
        ['test_zbl', 'ZBLPotTest.cpp', 'au_si'],
        ['test_extpot', 'ExtPotTest.cpp', 'extpot'],
        ['test_proj_rot_trans', 'ProjectOutRotTransTest.cpp', ''],
        ['test_epicenters', 'EpiCentersTest.cpp', 'Pt_Heptamer_FrozenLayers'],
        ['test_storage_order', 'StorageOrderTest.cpp', ''],
        ['test_helpers', 'HelperFunctionsTest.cpp', ''],
        ['test_lammps_loader', 'LammpsLoaderTest.cpp', ''],
        ['test_pot_registry', 'PotRegistryTest.cpp', ''],
        ['test_neb_lj', 'NEBMorseTest.cpp', 'neb_morse'],
        ['test_neb_regression', 'NEBRegressionTest.cpp', 'neb_morse'],
        ['test_dimer', 'DimerTest.cpp', 'neb_morse'],
        ['test_saddle', 'SaddleSearchTest.cpp', 'neb_morse'],
        ['test_confileio', 'ConFileIOTest.cpp', 'neb_morse'],
        ['test_minimization', 'MinimizationTest.cpp', 'neb_morse'],
        ['test_optimizer', 'OptimizerTest.cpp', ''],
        ['test_dynamics', 'DynamicsTest.cpp', 'neb_morse'],
        ['test_jobs', 'JobIntegrationTest.cpp', 'neb_morse'],
        ['test_hessian', 'HessianTest.cpp', 'neb_morse'],
        ['test_si_pot', 'SiPotTest.cpp', 'si_diamond'],
        ['test_eam_al', 'EAMAlTest.cpp', 'al_fcc'],
        ['test_fehe', 'FeHeTest.cpp', 'fe_bcc'],
        ['test_geom', 'GeometryAnalysisTest.cpp', 'neb_morse'],
        ['test_params_json', 'ParametersJSONTest.cpp', ''],
        ['test_bondboost', 'BondBoostTest.cpp', 'neb_morse'],
        ['test_emt_cu', 'EMTCuTest.cpp', 'cu_fcc'],
        ['test_neb_subsystem', 'NEBSubsystemTest.cpp', ''],
    ]
    if get_option('with_ase')
        test_array += [['test_ase_pot', 'ASEPotTest.cpp', 'ase_pot']]
    endif
    if get_option('with_metatomic')
        test_array += [['test_mta', 'MetatomicTest.cpp', 'lj38']]
    endif
    # XTB tests are always built; SKIP at runtime when libxtb isn't on
    # LD_LIBRARY_PATH (mirrors LAMMPS / ARTn / IRA tests).
    test_array += [['test_xtb', 'XTBTest.cpp', 'sulfolene']]
    test_array += [['test_cineb_xtb', 'CINEBXTBTest.cpp', 'cineb_xtb', 120]]
    if get_option('with_serve')
        test_array += [['test_serve_spec', 'ServeSpecParseTest.cpp', '']]
    endif
    # ARTn / IRA tests are always built; they SKIP at runtime when the
    # respective .so is not on LD_LIBRARY_PATH (see is_loaded() guards).
    test_array += [['test_artn', 'ARTnTest.cpp', 'neb_morse']]
    test_array += [['test_ira', 'IRATest.cpp', 'saddle_search']]
    foreach test : test_array
        test_timeout = test.get(3, 30)
        test(
            test.get(0),
            executable(
                test.get(0),
                sources: [
                    'unit_tests/' + test.get(1),
                    'thirdparty/catch2/catch_amalgamated.cpp',
                ],
                dependencies: test_deps,
                include_directories: _incdirs,
                cpp_args: test_args,
                link_with: _linkto,
                build_rpath: ':'.join(_runtime_rpath),
            ),
            workdir: meson.project_source_root() + '/client/unit_tests/data/systems/' + test.get(
                2,
            ),
            timeout: test_timeout,
        )
    endforeach
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
            link_with: _linkto,
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
        link_with: _linkto,
    )
    benchmark('bench_pot', bench_pot_exe,
        args: ['[.benchmark]'],
        workdir: meson.project_source_root() + '/client/unit_tests/data/systems/neb_morse',
        timeout: 120,
    )
endif
