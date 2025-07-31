potaids = library(
    'potaids',
    sources: files('ParseTOML.cc', 'PotHelpers.cc', 'PotentialCache.cc'),
    include_directories: _incdirs,
    dependencies: _deps,
    link_with: _linkto,
    cpp_args: _args,
    install: true,
)
_linkto += potaids

# subdir('EAM')
# subdir('EMT')
subdir('ExtPot')
# subdir('IMD')
subdir('LJ')
# subdir('LJCluster')
subdir('Morse')
# subdir('NewPot')

potentials = [
    # eam,
    # emt,
    extpot,
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

if get_option('with_vasp')
    if host_system != 'windows'
        subdir('VASP')
        potentials += [vasp]
    endif
endif

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

if get_option('with_lammps')
    subdir('LAMMPS')
    potentials += [lammps_pot]
endif

if py_embed
    _args += ['-DEMBED_PYTHON']
    # subdir('potentials/PyAMFF')
    # pyamff,
    # subdir('potentials/QSC')
    # potentials += [ qsc ]
    if get_option('with_ase')
        subdir('potentials/ASE')
        potentials += [ase]
        # TODO(rg): Annoying, use WITH_ style
        _args += ['-DASE_POT']
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
    subdir('Aluminum')
    subdir('EDIP')
    subdir('FeHe')
    subdir('Lenosky')
    subdir('SW')
    subdir('Tersoff')
    if get_option('with_water')
        subdir('Water')
        potentials += [water]
        subdir('Water_H')
        potentials += [water_h]
        subdir('Water_Pt')
        potentials += [water_pt]
    endif
    subdir('CuH2')
    potentials += [aluminum, cuh2, edip, fehe, lenosky, sw, tersoff]
    _args += ['-DWITH_FORTRAN']
endif

if get_option('with_ase_orca')
    # TODO: Cleanup, used for Catlearn too
    # Embedding the interpreter
    pyb11f_deps = [python_dep, dependency('pybind11')]
    _deps += [pyb11f_deps]
    subdir('ASE_ORCA')
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
    subdir('ASE_NWCHEM')
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
    LIB_TORCH_INCL_PATH_ATEN = LIB_TORCH_INCL_PATH / 'ATen'
    LIB_TORCH_INCL_PATH_C10 = LIB_TORCH_INCL_PATH / 'c10'
    LIB_TORCH_INCL_PATH_C10_UTIL = LIB_TORCH_INCL_PATH / 'c10' / 'util'

    LIB_TORCH_INCLUDE_LIST = [
        LIB_TORCH_INCL_PATH,
        LIB_TORCH_INCL_PATH_API,
        LIB_TORCH_INCL_PATH_ATEN,
        LIB_TORCH_INCL_PATH_C10,
        LIB_TORCH_INCL_PATH_C10_UTIL,
    ]


    lib_torch_list = ['c10', 'torch', 'torch_cpu', 'torch_global_deps']
    tdeps = []
    foreach lib_name : lib_torch_list
        tdeps += cppc.find_library(lib_name, dirs: [LIB_TORCH_LIB_PATH])
    endforeach

    torch_deps = declare_dependency(
        dependencies: tdeps,
        include_directories: include_directories(LIB_TORCH_INCLUDE_LIST),
    )

    if get_option('pip_metatomic')
        torch_ver = get_option('torch_version')
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
            include_directories: include_directories(
                py.get_install_dir() / 'metatensor' / 'include',
                py.get_install_dir() / 'metatensor' / 'torch' / f'torch-@torch_ver@' / 'include',
                py.get_install_dir() / 'metatomic' / 'torch' / f'torch-@torch_ver@' / 'include',
            ),
        )
        vesin_dep = declare_dependency(
            dependencies: [
                cppc.find_library(
                    'vesin',
                    dirs: [py.get_install_dir() / 'vesin' / 'lib'],
                ),
            ],
            include_directories: include_directories(
                py.get_install_dir() / 'vesin' / 'include',
            ),
        )
    else
        message(
            'Looking for existing installations of metatensor, metatensor-torch, metatomic-torch and vesin',
        )
        message(
            'If this is flaky, pass pip_metatomic=True with the torch version and let pip handle everything',
        )
        # XXX(rg): Still needs the includes, setup pkgconf and use dependency() with pkgconf instead
        # xref: https://github.com/metatensor/metatensor/issues/930
        vesin_dep = cppc.find_library('vesin')
        mts_dep = [
            cppc.find_library('metatensor'),
            cppc.find_library('metatensor_torch'),
            cppc.find_library('metatomic_torch'),
        ]
    endif
    _deps += [vesin_dep, mts_dep, torch_deps]
    subdir('Metatomic')
    potentials += [metatomic_pot]
    _args += ['-DWITH_METATOMIC']
endif

if get_option('with_vasp')
    if host_system != 'windows'
        subdir('VASP')
        potentials += [vasp]
        _args += ['-DWITH_VASP']
    endif
endif


_linkto += potentials
