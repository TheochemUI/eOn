A Windows MSVC build of pyeonclient no longer passes GCC-only flags such as `-Wno-unused-parameter`. GNU ld still wraps readcon-db with `--no-as-needed`, and macOS still links that library by path.
