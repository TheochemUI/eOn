potaids = library('potaids',
                sources : files(
                  [
                    'PotHelpers.cc',
                    'ParseTOML.cc',
                  ]
                ),
                include_directories: _incdirs,
                dependencies : _deps,
                link_with : _linkto,
                cpp_args : _args,
                install : true)
_linkto += potaids

# subdir('EAM')
# subdir('EMT')
# subdir('ExtPot')
# subdir('IMD')
subdir('LJ')
# subdir('LJCluster')
subdir('Morse')
# subdir('NewPot')

potentials = [
               # eam,
               # emt,
               # extpot,
               # imd,
               lennard_jones,
               # lennard_jones_cluster,
               morse,
               # newpot
             ]

# if get_option('with_gprd')
#   # subdir('GPRPotential')
#   # potentials += [ gprpot ]
# endif

# if get_option('with_gp_surrogate')
#   _args += ['-DWITH_GP_SURROGATE']
#   eonclib_sources += ['GPSurrogateJob.cpp']
#     if get_option('with_catlearn')
#       # TODO: Cleanup, used for ase_orca too
#       # Embedding the interpreter
#       pyb11f_deps = [ dependency('pybind11'),
#                       declare_dependency(link_args: '-lstdc++') ]
#       _deps += [ pyb11f_deps ]
#       subdir('CatLearnPot')
#       potentials += [ catlearnpot ]
#       _args += ['-DWITH_CATLEARN']
#     endif
# endif

# if get_option('with_vasp')
#   if host_system != 'windows'
#     subdir('VASP')
#     potentials += [ vasp ]
#     _args += ['-DWITH_VASP']
#   endif
# endif

# if get_option('with_water')
#   subdir('Water')
#   subdir('Water_Pt')
#   potentials += [ water, water_pt ]
#   _args += ['-DWITH_WATER']
# endif

# if get_option('with_ams')
#     # Optional, only for AMS
#     absl_dep = dependency('absl', method : 'cmake',
#       modules : [
#         'absl::base',
#         'absl::strings',
#         'absl::numeric',
#         'absl::utility',
#       ],
#       required : true,
#     )
#     boost_ams_dep = dependency('boost',
#       modules : [
#         'headers',
#         'asio',
#         'iostreams',
#       ],
#       required : true,
#     )
#     ams_dep = [ boost_ams_dep, absl_dep, fmt_dep  ]
#   _deps += ams_dep
#   # End AMS
#   subdir('AMS')
#   subdir('AMS_IO')
#   potentials += [ams, ams_io]
#   _args += ['-DWITH_AMS']
# endif

# if get_option('with_xtb')
#   xtb_proj = subproject('xtb')
#   xtb_dep = xtb_proj.get_variable('xtb_dep')
#   _deps += [xtb_dep]
#   subdir('XTBPot')
#   potentials += xtb_eon
#   _args += ['-DWITH_XTB']
# endif

# if get_option('with_lammps')
#   subdir('LAMMPS')
#   potentials += [ lammps_pot ]
#   _args += ['-DLAMMPS_POT']
# endif

# if get_option('with_python')
#   _args += ['-DWITH_PYTHON']
#   # subdir('PyAMFF')
#   # pyamff,
#   subdir('QSC')
#   potentials += [ qsc ]
#   if get_option('with_ase')
#     subdir('ASE')
#     potentials += [ ase ]
#     # TODO(rg): Annoying, use WITH_ style
#     _args += ['-DASE_POT']
#   endif
# endif

# if get_option('with_mpi')
#   subdir('MPIPot')
#   potentials += [ mpipot ]
#   _args += ['-DEONMPI']
#   mpi_dep = dependency('mpi')
#   _deps += [ mpi_dep ]
# endif

if get_option('with_fortran')
  subdir('Aluminum')
  subdir('EDIP')
#   subdir('FeHe')
#   subdir('Lenosky')
#   subdir('SW')
#   subdir('Tersoff')
#   if get_option('with_water')
#     subdir('Water_H')
#     potentials += [water_h]
#   endif
  subdir('CuH2')
  potentials += [
                 aluminum,
                 cuh2,
                 edip,
  ]
#   fehe, lenosky, sw, tersoff]
  _args += ['-DWITH_FORTRAN']
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
  subdir('ASE_ORCA')
  potentials += [ aseorca ]
  _args += ['-DWITH_ASE_ORCA']
endif

_linkto += potentials
