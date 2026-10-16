Linux keeps `--no-as-needed` around libreadcon_db. macOS records the dylib with `-needed_library`. Windows links the shared readcon-db library through the import library the subproject builds.
