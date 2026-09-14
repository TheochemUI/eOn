pyeonclient ConFrame export uses mkstemps (exclusive) and releases the GIL while writing. The guessable world-writable temp name is gone.
