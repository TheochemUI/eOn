# TODO(rg): Split this into pure python server and C++ client
# Most of this is only actually for the client
project('eon',
        'cpp',
        version : '2.1',
        default_options : [
          'warning_level=3',
          'cpp_std=c++17',
          # 'buildtype=debugoptimize', # debug=true, optimize=2
          'debug=true',
          'optimization=3',
        ])
# IMPORTANT!! warning_level=3 passes -fimplicit-none
# eonclient needs implicit typing!!


host_system = host_machine.system()

_args = [] # Extra arguments
_fargs = [] # Fortran arguments
_deps = [] # Dependencies
_linkto = [] # All the sub-libraries
_incdirs = [] # All the includes

subdir('client')

# eon main package
py = import('python').find_installation(pure: false)
py.install_sources([
  'eon/akmc.py',
  'eon/akmcstatelist.py',
  'eon/akmcstate.py',
  'eon/askmc.py',
  'eon/atoms.py',
  'eon/basinhopping.py',
  'eon/communicator.py',
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
  pure: false, # install next to compiled extension
  subdir: 'eon'
)

py.install_sources([
  'eon/migrator/config.py',
  'eon/migrator/config.yaml',
],
  pure: false,
  subdir: 'eon/migrator'
)
