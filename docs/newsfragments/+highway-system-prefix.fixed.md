The build no longer takes Highway from a system CMake package under `/usr`, which put `/usr/include` ahead of the compiler's own headers and broke conda builds. The wrap is used instead.
