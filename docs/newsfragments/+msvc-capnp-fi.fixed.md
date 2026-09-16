MSVC Cap'n Proto builds force-include a guard that undefines
windows.h `interface` after Meson's cl.exe sanity check.
