potss = ssmod.source_set()
potss.add(files('ParseTOML.cc', 'PotHelpers.cc'))

# subdir('EAM')
subdir('EMT')
subdir('ExtPot')
# subdir('IMD')
potss.add(files('LJ/LJ.cpp'))
# subdir('LJCluster')
potss.add(files('Morse/Morse.cpp'))
# subdir('NewPot')
potss.add(when: 'WITH_VASP', if_true: files('VASP/VASP.cpp'))
subdir('LAMMPS')
potss.add(
    when: 'WITH_FORTRAN',
    if_true: files(
# Aluminum needs help, since it has commonblks
# 'Aluminum/fofrhoDblexp.f',
# 'Aluminum/dfrhoDblexp.f',
# 'Aluminum/potinit.f',
# 'Aluminum/alforces.f',
# 'Aluminum/Aluminum.cpp',
# 'Aluminum/gagafeDblexp.f',
# 'Aluminum/embedenergy.f',
# 'Aluminum/sumembforce.f',
        'CuH2/CuH2.cpp',
        'CuH2/eam_dat.f90',
        'CuH2/eam_isoc.f90',
        'CuH2/eamroutines.f90',
        'EDIP/EDIP.cpp',
        'EDIP/edipFortran.f90',
        'FeHe/FeHe.cpp',
        'FeHe/feforce.f',
        'Lenosky/Lenosky.cpp',
        'Lenosky/lenoskyFortran.f90',
        'SW/SW.cpp',
        'SW/SWFortran.f90',
        'Tersoff/Tersoff.cpp',
        'Tersoff/tersoffFortran.f90',
    ),
)
potss.add(
    when: ['WITH_FORTRAN', 'WITH_WATER'],
    if_true: files(
        'Water/Water.cpp',
        'Water/ccl.cpp',
        'Water/potential_base.cpp',
        'Water/spce_ccl.cpp',
        'Water/tip4p_ccl.cpp',
        'Water_H/Tip4p_H.cpp',
        'Water_H/potH_H2O.f90',
        'Water_Pt/Tip4p_Pt.cpp',
        'Water_Pt/zhu_philpott.cpp',
        'Water_Pt/zhu_philpott_parameters.cpp',
    ),
)

# TODO(rg): Port to sourceset
# if get_option('with_gprd')
#   # subdir('GPRPotential')
#   # potentials += [ gprpot ]
# endif

# TODO(rg): Port to sourceset
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

if get_option('with_ase_orca')
# TODO: Cleanup, used for Catlearn too
    # Embedding the interpreter
    py = import('python').find_installation()
    python_dep = py.dependency(embed: true, required: true)
    pyb11f_deps = [
        python_dep,
        dependency('pybind11'),
        declare_dependency(link_args: '-lstdc++'),
    ]
    _deps += [pyb11f_deps]
    subdir('ASE_ORCA')
endif
