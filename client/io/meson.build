ecio_srcs = [
  'IOHelpers.cc',
  'IOBase.cc',
  'ConBase.cc',
  'ConWriter.cc',
  'ConvelWriter.cc',
  'XYZWriter.cc',
  'TibbleWriter.cc',
]

ecio = library(
  'ecio',
  sources : files(ecio_srcs),
  include_directories: _incdirs,
  dependencies : _deps,
  link_with : _linkto,
  cpp_args : _args,
  install : true
)
_linkto += ecio
