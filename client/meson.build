_incdirs += [ include_directories('thirdparty') ]

add_languages('c', required: true)
cc = meson.get_compiler('c')
cppc = meson.get_compiler('cpp')

use_fortran = false
if get_option('with_fortran') or get_option('with_xtb') or get_option('with_cuh2')
  use_fortran = true
  add_languages('fortran', required: true)
  fc = meson.get_compiler('fortran')
  _fargs += fc.get_supported_arguments([
    '-fno-implicit-none',
    ])
endif

if get_option('use_mkl')
    mkldep = dependency('mkl-dynamic-ilp64-iomp', required: true)
    add_project_arguments('-DEIGEN_USE_MKL_ALL', language: 'cpp')
    add_project_link_arguments('-lmkl_intel_lp64', '-lmkl_core', '-lmkl_sequential', language: 'cpp')
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

  incdir_numpy = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
    check : true
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

_args += cppc.get_supported_arguments([
  # C++ only
  '-Wno-class-memaccess',
  '-Wno-unused-function',
  # Common
  '-Wno-reorder',
  '-Wno-sign-compare',
  '-Wno-unused-but-set-parameter',
  '-Wno-unused-but-set-variable',
  '-Wno-unused-parameter',
  '-Wno-unused-result',
  '-Wno-unused-variable',
]
)

# TODO(rg): use fc.get_supported_arguments and add more compiler support
if use_fortran == true
    _fargs += fc.get_supported_arguments([
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
    ])
endif

# Debugging, should be offloaded from here eventually
_args += [
  '-g3',
  '-O0',
]

# Safer, saner default, can be turned off in release
if get_option('buildtype') != 'release'
  _args += [
    '-DEON_CHECKS=TRUE',
  ]
endif


# Get the version and build date
version_script = find_program('get_version.py')
version_output = run_command(version_script).stdout().strip()
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
  }
)

# Dependencies
if not is_windows
  # libm for Unix systems
  m_dep = cppc.find_library('m', required: false)
  _deps += m_dep
  # For building with clang
  _deps += [declare_dependency(link_args: '-lstdc++')]
endif

threads_dep = dependency('threads',
                         required: true)
_deps += threads_dep

eigen_dep = dependency('Eigen3', method : 'cmake', modules : ['Eigen3::Eigen'], required : true)
_deps += [ eigen_dep ]

fmt_dep = dependency('fmt',
                     required : true)
_deps += fmt_dep

spdlog_dep = dependency('spdlog', required : true)
_deps += spdlog_dep
_args += ['-DSPDLOG_FMT_EXTERNAL']

eoncbase_sources = [
  'INIFile.cpp',
]

eoncbase = library('eoncbase',
  sources : eoncbase_sources,
  dependencies : _deps,
  cpp_args : _args,
  include_directories: _incdirs,
  install : true)

_linkto += eoncbase

# Potentials
subdir('potentials/EAM')
subdir('potentials/EMT')
subdir('potentials/ExtPot')
subdir('potentials/IMD')
subdir('potentials/LJ')
subdir('potentials/LJCluster')
subdir('potentials/Morse')
subdir('potentials/NewPot')

potentials = [
               eam,
               emt,
               extpot,
               imd,
               lennard_jones,
               lennard_jones_cluster,
               morse,
               newpot
             ]

eonclib_sources = [
           'Parameters.cpp',
           'Optimizer.cpp',
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
           'TADJob.cpp',
           'ProcessSearchJob.cpp',
           'Bundling.cpp',
           'NudgedElasticBand.cpp',
           'MonteCarloJob.cpp',
           'DynamicsJob.cpp',
           'MonteCarlo.cpp',
           'Hessian.cpp',
           'NudgedElasticBandJob.cpp',
           'DynamicsSaddleSearch.cpp',
           'HelperFunctions.cpp',
           # 'StringHelpers.cc', # Template
           'MatrixHelpers.hpp', # Template
           'Dimer.cpp',
           'Dynamics.cpp',
           'GlobalOptimizationJob.cpp',
           'BiasedGradientSquaredDescent.cpp',
           'SafeHyperJob.cpp',
           'MinimizationJob.cpp',
           'Quickmin.cpp',
           'ParallelReplicaJob.cpp',
           'Potential.cpp',
           'SurrogatePotential.cpp', # Part of the interface
           'BasinHoppingJob.cpp',
           'FIRE.cpp',
           'EpiCenters.cpp',
           'SaddleSearchJob.cpp',
           'BasinHoppingSaddleSearch.cpp',
]

eonclient_sources = [
           'fpe_handler.cpp',
           'ClientEON.cpp',
           'CommandLine.cpp',
]

# ------------------- Options

if get_option('with_gprd')
    gprd_opts = {
        'default_library': 'static',
        'use_mkl': get_option('use_mkl'),
        'use_hdf5': true,
        'do_install': true,
        'buildtype': 'release',
    }
  libgprd_proj = subproject('gpr_optim',
                            default_options: gprd_opts)
  libgprd = libgprd_proj.get_variable('libgprd_dep')
  gprd_deps = libgprd_proj.get_variable('gprd_deps')
  # subdir('potentials/GPRPotential')
  # potentials += [ gprpot ]
  _args += ['-DWITH_GPRD']
  eonclib_sources += [ 'AtomicGPDimer.cpp', 'GPRHelpers.cpp' ]
  _deps += [ libgprd, gprd_deps ]
endif

if get_option('with_gp_surrogate')
  _args += ['-DWITH_GP_SURROGATE']
  eonclib_sources += ['GPSurrogateJob.cpp']
    if get_option('with_catlearn')
      # TODO: Cleanup, used for ase_orca too
      # Embedding the interpreter
      pyb11f_deps = [ dependency('pybind11') ]
      _deps += [ pyb11f_deps ]
      subdir('potentials/CatLearnPot')
      potentials += [ catlearnpot ]
      _args += ['-DWITH_CATLEARN']
    endif
endif

if get_option('with_vasp')
  if host_system != 'windows'
    subdir('potentials/VASP')
    potentials += [ vasp ]
    _args += ['-DWITH_VASP']
  endif
endif

if get_option('with_water')
  subdir('potentials/Water')
  subdir('potentials/Water_Pt')
  potentials += [ water, water_pt ]
  _args += ['-DWITH_WATER']
endif

if get_option('with_ams')
    # Optional, only for AMS
    absl_dep = dependency('absl', method : 'cmake',
      modules : [
        'absl::base',
        'absl::strings',
        'absl::numeric',
        'absl::utility',
      ],
      required : true,
    )
    boost_ams_dep = dependency('boost',
      modules : [
        'headers',
        'asio',
        'iostreams',
      ],
      required : true,
    )
    ams_dep = [ boost_ams_dep, absl_dep, fmt_dep  ]
  _deps += ams_dep
  # End AMS
  subdir('potentials/AMS')
  subdir('potentials/AMS_IO')
  potentials += [ams, ams_io]
  _args += ['-DWITH_AMS']
endif

if get_option('with_xtb')
  xtb_proj = subproject('xtb')
  xtb_dep = xtb_proj.get_variable('xtb_dep')
  _deps += [xtb_dep]
  subdir('potentials/XTBPot')
  potentials += xtb_eon
  _args += ['-DWITH_XTB']
endif

if get_option('with_lammps')
  subdir('potentials/LAMMPS')
  potentials += [ lammps_pot ]
  _args += ['-DLAMMPS_POT']
endif

if get_option('with_python')
  _args += ['-DWITH_PYTHON']
  # subdir('potentials/PyAMFF')
  # pyamff,
  # subdir('potentials/QSC')
  # potentials += [ qsc ]
  if get_option('with_ase')
    subdir('potentials/ASE')
    potentials += [ ase ]
    # TODO(rg): Annoying, use WITH_ style
    _args += ['-DASE_POT']
  endif
endif

if get_option('with_mpi')
  subdir('potentials/MPIPot')
  potentials += [ mpipot ]
  _args += ['-DEONMPI']
  mpi_dep = dependency('mpi')
  _deps += [ mpi_dep ]
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
  pyb11f_deps = [ python_dep,
                  dependency('pybind11') ]
  _deps += [ pyb11f_deps ]
  subdir('potentials/ASE_ORCA')
  potentials += [ aseorca ]
  _args += ['-DWITH_ASE_ORCA']
endif

if get_option('with_ase_nwchem')
  # TODO: Cleanup, used for Catlearn too
  # Embedding the interpreter
  py = import('python').find_installation()
  python_dep = py.dependency(embed: true, required: true)
  pyb11f_deps = [ python_dep,
                  dependency('pybind11'),
                  declare_dependency(link_args: '-lstdc++') ]
  _deps += [ pyb11f_deps ]
  subdir('potentials/ASE_NWCHEM')
  potentials += [ asenwchem ]
  _args += ['-DWITH_ASE_NWCHEM']
endif

_linkto += potentials

# --------------------- Library

eclib = library('eonclib',
                sources : eonclib_sources,
                include_directories: _incdirs,
                dependencies : _deps,
                link_with : _linkto,
                cpp_args : _args,
                link_args: _linkargs,
                install : true)

_linkto += eclib

# Export variable
eclib_dep = declare_dependency(
  link_with: [eclib, _linkto],
  include_directories: _incdirs,
  compile_args: _fargs + _args,
  dependencies: _deps,
)

# ---------------------- Executable

ec = executable('eonclient',
                sources : eonclient_sources,
                include_directories: _incdirs,
                dependencies : _deps,
                link_with : _linkto,
                cpp_args : [_args, '-w'],
                link_args: _linkargs,
                install : true)

if get_option('with_tests')
if get_option('with_python')
  # ec1 = executable('testPot',
  #                  'testpot.cpp',
  #                  include_directories: _incdirs,
  #                  dependencies : _deps,
  #                  link_with : _linkto,
  #                 )
endif

ec2 = executable('testWriteR',
                'testWriteR.cpp',
                include_directories: _incdirs,
                dependencies : _deps,
                link_with : _linkto,
                )
endif
# ------------------------ Tests

if get_option('with_tests')
    test_args = _args
    test_deps = _deps
    testMain = library(
      'TestMain',
      ['gtests/TestMain.cpp', 'thirdparty/catch2/catch_amalgamated.cpp'],
      dependencies : test_deps,
      include_directories: _incdirs,
      cpp_args: test_args,
      link_with: _linkto,
    )
    test_array = [#
      # ['test_impldim', 'ImpDimerTest.cpp', '/gtests/data/saddle_search'],
      # ['strparse_run', 'StringHelpersTest.cpp', ''],
      ['test_matter', 'MatterTest.cpp', 'sulfolene'],
      ['test_pot', 'PotTest.cpp', 'sulfolene'],
      ['approve_params', 'ApproveParams.cpp', ''],
      # ['neb_test', 'NEBTest.cpp', '/gtests/data/saddle_search'],
    ]
    foreach test : test_array
      test(test.get(0),
           executable(test.get(0),
              sources : ['gtests/'+test.get(1), 'thirdparty/catch2/catch_amalgamated.cpp'],
              dependencies : test_deps,
              include_directories: _incdirs,
              cpp_args: test_args,
              link_with: _linkto + testMain,
                     ),
            workdir : meson.source_root() + '/client/gtests/data/systems/' + test.get(2)
          )
    endforeach
endif
