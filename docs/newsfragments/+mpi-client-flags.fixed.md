An MPI build of `eonclient` reads its command-line flags (`--version`, `--help`, one-shot jobs) unless an eOn server launched the rank. It used to ignore every flag and run as a client.
