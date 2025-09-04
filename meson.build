# TODO(rg): Split this into pure python server and C++ client
# Most of this is only actually for the client
project(
    'eon',
    'cpp',
    version: run_command(['tools/gitversion.py'], check: true).stdout().strip(),
    meson_version: '>= 1.6.0',
    default_options: [  #
        'buildtype=debugoptimized',
        'warning_level=0',
        'cpp_std=c++17',
    ],
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
subdir('eon')
