project('eonclient', 'cpp',
  version : '0.1',
  default_options : ['warning_level=2', 'cpp_std=c++17'])
# IMPORTANT!! warning_level=3 passes -fimplicit-none
# EON needs implicit typing!!

host_system = host_machine.system()
eon_extra_args = []

add_languages('c', required: true)
cc = meson.get_compiler('c')
cppc = meson.get_compiler('cpp')

# Add conditionals
if host_system == 'darwin'
  eon_extra_args += ['-DOSX=TRUE']
  # Workaround for
  # ../MinModeSaddleSearch.cpp:173:25: error: aligned allocation function of
  # type 'void *(unsigned long, enum std::align_val_t)' is only available on macOS 10.13 or newer
  add_global_arguments(['-faligned-allocation'], language: 'cpp')
endif

# Dependencies
eigen_dep = dependency('Eigen3', method : 'cmake', modules : ['Eigen3::Eigen'], required : true)
eon_deps = [ eigen_dep ]

fmt_dep = dependency('fmt', method : 'cmake',  modules : ['fmt::fmt'],  required : true)
eon_deps += fmt_dep

# Testing only
gtest_dep = dependency('gtest', main : true, required : true)
gmock_dep = dependency('gmock', required : false)

# Optional, only for AMS
absl_dep = dependency('absl', method : 'cmake',
  modules : [
    'absl::base',
    'absl::strings',
    'absl::numeric',
    'absl::utility',
  ],
  required : false,
)

boost_ams_dep = dependency('boost',
  modules : [
    'headers',
    'asio',
    'iostreams',
  ],
  required : false,
)

ams_dep = [ boost_ams_dep, absl_dep, fmt_dep  ]

# End AMS

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

eon_pot = declare_dependency(link_with : potentials)

eonclib_sources = [
           'Optimizer.cpp',
           'PrefactorJob.cpp',
           'LBFGS.cpp',
           'Parameters.cpp',
           'ReplicaExchangeJob.cpp',
           'BondBoost.cpp',
           'Job.cpp',
           'INIFile.cpp',
           'GlobalOptimization.cpp',
           'Log.cpp',
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
           'StringHelpers.cc', # Template
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
           'BaseStructures.cpp',
]

eonclient_sources = [
           'ClientEON.cpp',
           'CommandLine.cpp',
]

# ------------------- Options

if get_option('with_gprd')
  libgprd_proj = subproject('gprdimer')
  gprd_dep = libgprd_proj.get_variable('libgprd_dep')
  subdir('potentials/GPRPotential')
  potentials += [ gprpot ]
  eon_extra_args += ['-DWITH_GPRD']
  eonclient_sources += [ 'AtomicGPDimer.cpp', 'GPRHelpers.cpp' ]
  eon_deps += [ gprd_dep ]
endif

if get_option('with_vasp')
  if host_system != 'windows'
    subdir('potentials/VASP')
    potentials += [ vasp ]
    eon_extra_args += ['-DWITH_VASP']
  endif
endif

if get_option('with_water')
  subdir('potentials/Water')
  subdir('potentials/Water_Pt')
  potentials += [ water, water_pt ]
  eon_extra_args += ['-DWITH_WATER']
endif

if get_option('with_ams')
  subdir('potentials/AMS')
  subdir('potentials/AMS_IO')
  potentials += [ams, ams_io]
  eon_extra_args += ['-DWITH_AMS']
endif

if get_option('with_lammps')
  subdir('potentials/LAMMPS')
  potentials += [ lammps_pot ]
  eon_extra_args += ['-DLAMMPS_POT']
endif

if get_option('with_python')
  subdir('potentials/PyAMFF')
  subdir('potentials/QSC')
  potentials += [ pyamff, qsc ]
  eon_extra_args += ['-DWITH_PYTHON']
endif

if get_option('with_mpi')
  subdir('potentials/MPIPot')
  potentials += [ mpipot ]
  eon_extra_args += ['-DEONMPI']
  mpi_dep = dependency('mpi')
  eon_deps += [ mpi_dep ]
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
  eon_extra_args += ['-DWITH_FORTRAN']
endif

if get_option('with_cuh2')
  # Implicitly needs Fortran
  add_languages('fortran', required: true)
  fc = meson.get_compiler('fortran')
  subdir('potentials/CuH2')
  potentials += [cuh2]
  eon_extra_args += ['-DCUH2_POT']
endif

# --------------------- Library

eclib = library('eonclib',
                sources : eonclib_sources,
                dependencies : eon_deps,
                link_with : potentials,
                cpp_args : eon_extra_args,
                install : true)

# ---------------------- Bindings

if get_option('with_pybindings')
  # Implicitly needs Fortran
  py = import('python').find_installation()
  python_dep = py.dependency()
  pyb11f_deps = [ python_dep ]


  pyb11f_deps += [eon_deps,
                  dependency('python3'),
                  dependency('pybind11'),
                  declare_dependency(link_args: '-lstdc++') ]

  pyb11f = py.extension_module('pyeonclient',
                               sources : [
                                 'pybinds/py_wrapper.cc',
                                 'pybinds/py_parameters.cc',
                               ],
                               dependencies : pyb11f_deps,
                               link_with : eclib,
                               cpp_args : eon_extra_args,
                               install : true)
endif

# ---------------------- Executable

ec = executable('eonclient',
                sources : eonclient_sources,
                dependencies : eon_deps,
                link_with : eclib,
                cpp_args : eon_extra_args,
                install : true)

ec1 = executable('testPot',
                'testpot.cpp',
                dependencies : eon_deps,
                link_with : eclib,
                )

ec2 = executable('testWriteR',
                'testWriteR.cpp',
                dependencies : eon_deps,
                link_with : eclib,
                )

# ------------------------ Tests

if get_option('with_tests')
  eon_extra_args += ['-DEONTEST'] # Unused
  eon_deps += [ gtest_dep, gmock_dep ]
test_array = [#
  ['Improved Dimer', 'impldim_run', 'ImpDimerTest.cpp', '/gtests/data/saddle_search'],
  ['String parser helpers', 'strparse_run', 'StringHelpersTest.cpp', ''],
  ['Matter converter', 'matter_run', 'MatterTest.cpp', '/gtests/data/systems/sulfolene'],
  ['Potential tests', 'test_pot', 'PotTest.cpp', '/gtests/data/systems/sulfolene'],
  ['NEB tests', 'neb_test', 'NEBTest.cpp', '/gtests/data/saddle_search'],
             ]
foreach test : test_array
  test(test.get(0),
       executable(test.get(1),
          sources : ['gtests/'+test.get(2)],
          dependencies : [ eon_deps, gtest_dep, gmock_dep ],
          link_with : eclib,
          cpp_args : eon_extra_args
                 ),
        workdir : meson.source_root() + test.get(3)
      )
endforeach
endif
