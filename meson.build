# TODO(rg): Split this into pure python server and C++ client
# Most of this is only actually for the client
project(
    'eon',
    'cpp',
    version: '3.0',
    default_options: ['warning_level=0', 'cpp_std=c++17'],
)
# IMPORTANT!! warning_level=3 passes -fimplicit-none
# eonclient needs implicit typing!!


host_system = host_machine.system()

_args = []  # Extra arguments
_fargs = []  # Fortran arguments
_deps = []  # Dependencies
_linkto = []  # All the sub-libraries
_linkargs = []  # Linker arguments
_incdirs = [include_directories('.')]  # All the includes

# These are actually from the client, but the py interpreter is determined here
py_modules = []
py_embed = false
need_ase = get_option('with_ase_orca') or get_option('with_ase_nwchem') or get_option(
    'with_ase',
) or get_option(
    'with_catlearn',
)
if need_ase
    py_modules += ['numpy', 'ase']
    py_embed = true
endif

if get_option('with_catlearn')
    py_modules += ['catlearn']
    py_embed = true
endif

if get_option('with_qsc')
    py_modules += ['cython']
    py_embed = true
endif

if get_option('pip_metatomic')
    py_modules += ['vesin', 'metatensor', 'metatensor.torch', 'metatomic.torch']
    # We don't need the runtime, so we don't change py_embed
endif

py = import('python').find_installation(pure: false, modules: py_modules)
if py_embed
    python_dep = py.dependency(embed: py_embed)
    _deps += python_dep
    # Safe to say that we need numpy as a dep if we are embedding
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

# Client binary
subdir('client')

# eon main package
py.install_sources(
    [
        'eon/akmc.py',
        'eon/akmcstatelist.py',
        'eon/akmcstate.py',
        'eon/askmc.py',
        'eon/atoms.py',
        'eon/basinhopping.py',
        'eon/communicator.py',
        'eon/config.py',
        'eon/config.yaml',
        'eon/displace.py',
        'eon/eon_kdb.py',
        'eon/escaperate.py',
        'eon/explorer.py',
        'eon/fileio.py',
        'eon/__init__.py',
        'eon/locking.py',
        'eon/server.py',
        'eon/schema.py',
        # TODO(rg): Spin mcamc out and use as a dependency
        # This will only use the python variant for now
        'eon/mcamc/__init__.py',
        'eon/mcamc/mcamc.py',
        'eon/movie.py',
        'eon/mpiwait.py',
        'eon/parallelreplica.py',
        'eon/prstatelist.py',
        'eon/prstate.py',
        'eon/recycling.py',
        'eon/statelist.py',
        'eon/state.py',
        'eon/superbasin.py',
        'eon/superbasinscheme.py',
        'eon/version.py',
    ],
    pure: false,  # install next to compiled extension
    subdir: 'eon',
)
