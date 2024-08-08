_incdirs += [include_directories('thirdparty')]

add_languages('c', required: true)
cc = meson.get_compiler('c')
cppc = meson.get_compiler('cpp')

use_fortran = false
if get_option('with_fortran') or get_option('with_xtb') or get_option(
    'with_cuh2',
)
    use_fortran = true
    add_languages('fortran', required: true)
    fc = meson.get_compiler('fortran')
endif

if get_option('with_python')
  # TODO(rg): Cleanup a bit
  # Embed is needed for surrogate, numpy by qsc
  if get_option('with_ase_orca')
    py = import('python').find_installation('python', modules: ['numpy', 'ase'])
  else
    py = import('python').find_installation('python', modules: ['numpy'])
  endif
  python_dep = py.dependency(embed: true)
  _deps += python_dep

    incdir_numpy = run_command(
        py,
        [
            '-c',
            'import os; os.chdir(".."); import numpy; print(numpy.get_include())',
        ],
        check: true,
    ).stdout().strip()
    inc_np = include_directories(incdir_numpy)
    _deps += declare_dependency(include_directories: inc_np)
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
    _args += ['-DIS_WINDOWS']
# Workaround used in INIFile.cpp
# ../11.4.0/include/c++/bits/basic_string.h:1059: std::__cxx11::basic_string<_CharT, _Traits, _Alloc>::reference std::__cxx11::basic_string<_CharT, _Traits, _Alloc>::operator[](std::__cxx11::basic_string<_CharT, _Traits, _Alloc>::size_type) [with _CharT = char; _Traits = std::char_traits<char>; _Alloc = std::allocator<char>; std::__cxx11::basic_string<_CharT, _Traits, _Alloc>::reference = char&; std::__cxx11::basic_string<_CharT, _Traits, _Alloc>::size_type = long unsigned int]: Assertion '__pos <= size()' failed.
endif

_args += cppc.get_supported_arguments(
    [
# C++ only
# '-Wno-class-memaccess',
# '-Wno-unused-function',
# Common
        '-Wno-unused-parameter',
        '-Wno-unused-but-set-parameter',
        '-Wfatal-errors',
        # From Eigen
        '-Wno-maybe-uninitialized',
        # From cachelot
        '-Wno-ignored-qualifiers',
    ],
)

# TODO(rg): use fc.get_supported_arguments and add more compiler support
if use_fortran == true
    _fargs += fc.get_supported_arguments(
        [
# Ported from C++
            '-Wno-unused-but-set-variable',
            '-Wno-unused-variable',
            # Absolutely necessary
            '-fno-implicit-none',
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


# Get the version and build date
version_script = find_program('get_version.py')
version_output = run_command(version_script, check: true).stdout().strip()
version_bits = version_output.split(',')
version = version_bits[0]
build_date = version_bits[1]
architecture = host_machine.cpu_family()
message(host_system)
message(architecture)

# Configure the header file
configure_file(
    input: 'version.h.in',
    output: 'version.h',
    configuration: {
        'VERSION': version,
        'BUILD_DATE': build_date,
        'OS_INFO': host_system,
        'ARCH': architecture,
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

eigen_dep = dependency(
    'Eigen3',
    method: 'cmake',
    modules: ['Eigen3::Eigen'],
    required: true,
)
_deps += [eigen_dep]

fmt_dep = dependency('fmt', required: true)
_deps += fmt_dep

spdlog_dep = dependency('spdlog', required: true)
_deps += spdlog_dep
_args += ['-DSPDLOG_FMT_EXTERNAL']

# Handle client flags
if get_option('with_lammps')
    _args += ['-DLAMMPS_POT']
endif

if get_option('with_fortran')
    _args += ['-DWITH_FORTRAN']
    if get_option('with_water')
        _args += ['-DWITH_WATER']
    endif
endif

if get_option('with_ase_orca')
    _args += ['-DWITH_ASE_ORCA']
endif

if get_option('with_vasp')
    if host_system != 'windows'
        _args += ['-DWITH_VASP']
    endif
endif

# Subprojects
cmake = import('cmake')
optv_cachelot = cmake.subproject_options()
optv_cachelot.add_cmake_defines(
    {
        'CMAKE_BUILD_TYPE': 'Release',
        'Boost_USE_STATIC_LIBS': false,
        'NO_SERVER': true,
        'CMAKE_POSITION_INDEPENDENT_CODE': true,
        'CMAKE_CXX_FLAGS': '-w',
        'CMAKE_C_FLAGS': '-w',
    },
)
cmcachelot = cmake.subproject('cachelot', options: optv_cachelot)
_deps += cmcachelot.dependency('cachelot')

optv_xxhashpp = cmake.subproject_options()
optv_xxhashpp.add_cmake_defines(
    {'CMAKE_BUILD_TYPE': 'Release', 'CMAKE_POSITION_INDEPENDENT_CODE': true},
)
cmxxhashpp = cmake.subproject('xxhashpp', options: optv_xxhashpp)
_deps += cmxxhashpp.dependency('xxhash_cpp')

# Potentials
subdir('potentials')

tomlib = library('toml', 'toml.cpp', cpp_args: ['-w'], install: true)
_deps += declare_dependency(link_with: tomlib)

ecparser_srcs = ['Parser.cc', 'EnvHelpers.cc']

ecparser = library(
    'ecparser',
    sources: files(ecparser_srcs),
    include_directories: _incdirs,
    dependencies: _deps,
    link_with: _linkto,
    cpp_args: _args,
    install: true,
)

_linkto += ecparser

_linkto += library(
    'echelpers',
    sources: files('HelperFunctions.cpp'),
    include_directories: _incdirs,
    dependencies: _deps,
    link_with: _linkto,
    cpp_args: _args,
    install: true,
)

subdir('matter')

# IO operations
subdir('io')

eonclib_sources = [
# 'Optimizer.cpp',
# 'PrefactorJob.cpp',
# 'LBFGS.cpp',
# 'ReplicaExchangeJob.cpp',
# 'BondBoost.cpp',
    'Parameters.cpp',  # TODO(rg) remove
    'Job.cpp',
    # 'GlobalOptimization.cpp',
    # 'LowestEigenmode.cpp',
    # 'MinModeSaddleSearch.cpp',
    # 'StructureComparisonJob.cpp',
    # 'SteepestDescent.cpp',
    # 'ImprovedDimer.cpp',
    'PointJob.cpp',
    # 'Prefactor.cpp',
    # 'ConjugateGradients.cpp',
    # 'FiniteDifferenceJob.cpp',
    # 'Lanczos.cpp',
    # 'HessianJob.cpp',
    # 'TADJob.cpp',
    # 'ProcessSearchJob.cpp',
    'Bundling.cpp',
    # 'NudgedElasticBand.cpp',
    # 'MonteCarloJob.cpp',
    # 'DynamicsJob.cpp',
    # 'MonteCarlo.cpp',
    # 'Hessian.cpp',
    # 'NudgedElasticBandJob.cpp',
    # 'DynamicsSaddleSearch.cpp',
    # # 'StringHelpers.cc', # Template
    # 'MatrixHelpers.hpp', # Template
    # 'Dimer.cpp',
    # 'Dynamics.cpp',
    # 'GlobalOptimizationJob.cpp',
    # 'BiasedGradientSquaredDescent.cpp',
    # 'SafeHyperJob.cpp',
    # 'MinimizationJob.cpp',
    # 'Quickmin.cpp',
    # 'ParallelReplicaJob.cpp',
    'Potential.cpp',
    # 'SurrogatePotential.cpp', # Part of the interface
    # 'BasinHoppingJob.cpp',
    # 'FIRE.cpp',
    # 'EpiCenters.cpp',
    # 'SaddleSearchJob.cpp',
    # 'BasinHoppingSaddleSearch.cpp',
]

eonclient_sources = ['fpe_handler.cpp', 'ClientEON.cpp', 'CommandLine.cpp']

# ------------------- Options

# if get_option('with_gprd')
#   libgprd_proj = subproject('gpr_optim')
#   gprd_dep = libgprd_proj.get_variable('libgprd_dep')
#   _args += ['-DWITH_GPRD']
#   eonclib_sources += [ 'AtomicGPDimer.cpp', 'GPRHelpers.cpp' ]
#   _deps += [ gprd_dep ]
# endif

# --------------------- Library

eclib = library(
    'eonclib',
    sources: eonclib_sources,
    include_directories: _incdirs,
    dependencies: _deps,
    link_with: _linkto,
    cpp_args: _args,
    link_args: _linkargs,
    install: true,
)

_linkto += eclib

# Export variable
eclib_dep = declare_dependency(
    link_with: [eclib, _linkto],
    include_directories: _incdirs,
    compile_args: _fargs + _args,
    dependencies: _deps,
)

# ---------------------- Executable

ec = executable(
    'eonclient',
    sources: eonclient_sources,
    include_directories: _incdirs,
    dependencies: _deps,
    link_with: _linkto,
    cpp_args: [_args],
    link_args: _linkargs,
    install: true,
)

if get_option('with_tests')
    if get_option('with_python')
# ec1 = executable('testPot',
        #                  'testpot.cpp',
        #                  include_directories: _incdirs,
        #                  dependencies : _deps,
        #                  link_with : _linkto,
        #                 )
    endif

# ec2 = executable('testWriteR',
#                 'testWriteR.cpp',
#                 include_directories: _incdirs,
#                 dependencies : _deps,
#                 link_with : _linkto,
#                 )
endif
# ------------------------ Tests

if get_option('with_tests')
    catch_cpp = files('thirdparty/catch2/catch_amalgamated.cpp')
    subdir('gtests')
endif
