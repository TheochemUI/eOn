Same-size ``Matter::resize`` keeps .con column-5 atom ids. A leftover
duplicate resize body was resetting them to ``0..N-1``.
