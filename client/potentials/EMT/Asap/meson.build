asap_files = [
  'NeighborList.cpp',
  'GhostAtoms.cpp',
  'CellList.cpp',
  'Vec.cpp',
  'Atoms.cpp',
  'Exception.cpp',
  'EMTRasmussenParameterProvider.cpp',
  'EMTDefaultParameterProvider.cpp',
  'EMT.cpp',
  'SuperCell.cpp',
]
# Don't want to deal with warnings from embedded lib
asap_args = [_args, '-w']
if get_option('with_asap_stacktrace')
  asap_files += ['stacktrace.c']
  asap_args = ['-DSTACKTRACE']
endif
asap = shared_library('asap', asap_files,
                      dependencies : _deps,
                      link_with : _linkto,
                      cpp_args : asap_args,
                      install : true)
