On Windows, ExtPot launches a suffix-less ext_pot wrapper with python.
cmd.exe does not honor a shebang, so a quoted path to the script was
not a runnable command.
