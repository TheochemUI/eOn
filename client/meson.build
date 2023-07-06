project('eonclient', 'cpp',
  version : '0.1',
  default_options : ['warning_level=2', 'cpp_std=c++17'])
# IMPORTANT!! warning_level=3 passes -fimplicit-none
# EON needs implicit typing!!


host_system = host_machine.system()

_args = [] # Extra arguments
_deps = [] # Dependencies
_linkto = [] # All the sub-libraries
_incdirs = [] # All the includes

# Someday these should be removed and things should be fixed..
_args += ['-Wno-unused-result',
         '-Wunused-but-set-parameter']

add_languages('c', required: true)
cc = meson.get_compiler('c')
cppc = meson.get_compiler('cpp')

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
                     version: ['<10'],
                     required : true)
_deps += fmt_dep

spdlog_dep = dependency('spdlog', required : true)
_deps += spdlog_dep
_args += ['-DSPDLOG_FMT_EXTERNAL']

eoncbase_sources = [
  'BaseStructures.cpp',
  'INIFile.cpp',
]

eoncbase = library('eoncbase',
  sources : eoncbase_sources,
  dependencies : _deps,
  cpp_args : _args,
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
           'BasinHoppingJob.cpp',
           'FIRE.cpp',
           'EpiCenters.cpp',
           'SaddleSearchJob.cpp',
           'BasinHoppingSaddleSearch.cpp',
]

eonclient_sources = [
           'ClientEON.cpp',
           'CommandLine.cpp',
]

# ------------------- Options

if get_option('with_gpr_optim')
  libgprd_proj = subproject('gpr_optim')
  gprd_dep = libgprd_proj.get_variable('libgprd_dep')
  _deps += [ gprd_dep ]
  _args += ['-DWITH_GPR_OPTIM']
  subdir('potentials/GPRPotential')
  potentials += [ gprpot ]
  eonclib_sources += [ 'AtomicGPDimer.cpp', 'GPRHelpers.cpp' ]
endif

if get_option('with_gp_surrogate')
  _args += ['-DWITH_GP_SURROGATE']
  eonclib_sources += ['SurrogatePotential.cpp', 'GPSurrogateJob.cpp']
    if get_option('with_catlearn')
      # TODO: Cleanup, used for ase_orca too
      # Embedding the interpreter
      py = import('python').find_installation()
      python_dep = py.dependency(embed: true)
      pyb11f_deps = [ python_dep,
                      dependency('pybind11'),
                      declare_dependency(link_args: '-lstdc++') ]
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
  # TODO: Clean this up and integrate with with_fortran
  add_languages('fortran', required: true)
  fc = meson.get_compiler('fortran')
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
  subdir('potentials/PyAMFF')
  subdir('potentials/QSC')
  potentials += [ pyamff, qsc ]
  _args += ['-DWITH_PYTHON']
endif

if get_option('with_mpi')
  subdir('potentials/MPIPot')
  potentials += [ mpipot ]
  _args += ['-DEONMPI']
  mpi_dep = dependency('mpi')
  _deps += [ mpi_dep ]
endif

if get_option('with_fortran')
  add_languages('fortran', required: true)
  fc = meson.get_compiler('fortran')
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
  # Implicitly needs Fortran
  add_languages('fortran', required: true)
  fc = meson.get_compiler('fortran')
  subdir('potentials/CuH2')
  potentials += [cuh2]
  _args += ['-DCUH2_POT']
endif

if get_option('with_ase_orca')
  # TODO: Cleanup, used for Catlearn too
  # Embedding the interpreter
  py = import('python').find_installation()
  python_dep = py.dependency(embed: true, required: true)
  pyb11f_deps = [ python_dep,
                  dependency('pybind11'),
                  declare_dependency(link_args: '-lstdc++') ]
  _deps += [ pyb11f_deps ]
  subdir('potentials/ASE_ORCA')
  potentials += [ aseorca ]
  _args += ['-DWITH_ASE_ORCA']
endif

_linkto += potentials

# --------------------- Library

eclib = library('eonclib',
                sources : eonclib_sources,
                dependencies : _deps,
                link_with : _linkto,
                cpp_args : _args,
                install : true)

_linkto += eclib

# ---------------------- Executable

ec = executable('eonclient',
                sources : eonclient_sources,
                dependencies : _deps,
                link_with : _linkto,
                cpp_args : _args,
                install : true)

ec1 = executable('testPot',
                'testpot.cpp',
                dependencies : _deps,
                link_with : _linkto,
                )

ec2 = executable('testWriteR',
                'testWriteR.cpp',
                dependencies : _deps,
                link_with : _linkto,
                )

# ------------------------ Tests

if get_option('with_tests')
  gtest_dep = dependency('gtest', main : true, required : true)
  gmock_dep = dependency('gmock', required : false)
  _args += ['-DEONTEST'] # Unused
  _deps += [ gtest_dep, gmock_dep ]
test_array = [#
  # ['Improved Dimer', 'test_impldim', 'ImpDimerTest.cpp', '/gtests/data/saddle_search'],
  # ['String parser helpers', 'strparse_run', 'StringHelpersTest.cpp', ''],
  ['Matter converter', 'test_matter', 'MatterTest.cpp', '/gtests/data/systems/sulfolene'],
  ['Potential tests', 'test_pot', 'PotTest.cpp', '/gtests/data/systems/sulfolene'],
  # ['NEB tests', 'neb_test', 'NEBTest.cpp', '/gtests/data/saddle_search'],
             ]
foreach test : test_array
  test(test.get(0),
       executable(test.get(1),
          sources : ['gtests/'+test.get(2)],
          dependencies : _deps,
          link_with : _linkto,
          cpp_args : _args
                 ),
        workdir : meson.source_root() + test.get(3)
      )
endforeach
endif
