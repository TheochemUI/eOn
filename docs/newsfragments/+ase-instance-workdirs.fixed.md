ASE-NWChem and ASE-ORCA calculators each get a private work directory instead
of `directory='.'`, so concurrent LocalInProcess jobs do not clobber each
other's scratch files. Calculator errors throw rather than abort the process.
