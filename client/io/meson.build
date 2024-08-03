ecio = library('ecio',
                sources : files(
                  [
                    'IOBase.cc',
                    'ConWriter.cc',
                  ]
                ),
                include_directories: _incdirs,
                dependencies : _deps,
                link_with : _linkto,
                cpp_args : _args,
                install : true)
_linkto += ecio
