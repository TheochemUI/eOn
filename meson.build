# TODO(rg): Split this into pure python server and C++ client
# Most of this is only actually for the client
project(
    'eon',
    'cpp',
    version: run_command(['tools/gitversion.py'], check: true).stdout().strip(),
    meson_version: '>= 1.8.0',
    default_options: [  #
        'buildtype=debugoptimized',
        'warning_level=0',
        'cpp_std=c++20',
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
# Project root (subprojects/...) + package include tree include/eon/
_incdirs = [include_directories('.'), include_directories('include')]

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

if get_option('pip_metatomic')
    py_modules += ['vesin', 'metatensor', 'metatensor.torch', 'metatomic.torch']
    # We don't need the runtime, so we don't change py_embed
endif

# On Windows, conda/pixi often put a PE binary at `python` (no .exe). Meson then
# tries to read it as a shebang script, hits UnicodeDecodeError on the MZ header
# ("Unusable script"), and find_installation fails the pip_metatomic modules.
# Prefer python.exe so ExternalProgram treats it as a real Windows binary.
_py_cmd = host_system == 'windows' ? 'python.exe' : 'python3'
# find_installation(..., modules:) only reports "missing modules" with no stderr.
# On Windows pip torch, imports need DLL directories registered; use the helper
# for a single stack check. Elsewhere check each module so ImportError is logged.
py = import('python').find_installation(_py_cmd, pure: false)
if py_modules.length() > 0
    if host_system == 'windows' and get_option('pip_metatomic')
        _r = run_command(
            py,
            meson.project_source_root() / 'tools' / 'win_pip_torch_path.py',
            '--import-check',
            check: false,
        )
        if _r.returncode() != 0
            error(
                'pip_metatomic Python stack not importable with '
                + py.full_path()
                + ':\n'
                + _r.stderr()
                + _r.stdout(),
            )
        endif
    else
        foreach _mod : py_modules
            _r = run_command(py, '-c', 'import ' + _mod, check: false)
            if _r.returncode() != 0
                error(
                    'Python module "' + _mod + '" is not importable with '
                    + py.full_path()
                    + ':\n'
                    + _r.stderr()
                    + _r.stdout(),
                )
            endif
        endforeach
    endif
endif
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
    _deps += declare_dependency(compile_args: '-I' + incdir_numpy)
endif

# Client binary / libs (+ optional pyeonclient)
subdir('client')

# Public headers: $prefix/include/eon/... (numpy/fmt package layout).
# Sources stay under client/; consumers use -I$prefix/include and
# #include "eon/Potential.h".
install_subdir(
    'include/eon',
    install_dir: get_option('includedir'),
)

# eon main package (AKMC server). Skip for standalone pyeonclient wheels.
if get_option('install_eon_server')
    subdir('eon')
endif
