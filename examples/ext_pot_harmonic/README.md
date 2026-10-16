# External harmonic potential

This standalone example runs eOn against a Python file-interface wrapper.
It needs no electronic-structure program, MPI runtime, model or
pseudopotentials. The potential is
`E = k/2 sum_i |r_i|^2`, with `k = 1 eV/Angstrom^2`, and the force is
`F_i = -k r_i`. It is an interface demonstration with a fixed origin.

Install [Pixi](https://pixi.sh) and run from this directory in a Linux
terminal (including Ubuntu under WSL):

```bash
curl -fsSL https://pixi.sh/install.sh | sh
export PATH="$HOME/.pixi/bin:$PATH"
pixi install
bash run.sh
```

The manifest installs the packaged eOn client and Python from conda-forge;
there is no source compilation. Keep the generated `pixi.lock` alongside
the manifest when sharing an exact environment, and recreate that
environment with `pixi install --locked`.

The supplied `pos.con` contains one silicon atom at `(1, 0, 0)` Angstrom
in a 20 Angstrom cubic box. `run.sh` creates a distinct directory under
`runs/`, copies the three inputs there, runs a point job and checks the
result. An optional first argument selects another work directory:

```bash
bash run.sh /absolute/path/to/scratch/harmonic-point
```

Expected values are **0.5 eV** energy and **(-1, 0, 0) eV/Angstrom** force.
The result check verifies energy and maximum force (1 eV/Angstrom), then
prints the record directory. `results.dat` and `eon.log` remain there.
An existing result is preserved and requires a different work directory.

On an HPC system, run `bash run.sh` inside a compute allocation. One CPU
is sufficient for this example. The files also run on Windows through
WSL; `run.sh` requires a shell and is not a PowerShell script.

`harmonic_extpot.py` reads `from_eon_to_extpot`: three box rows in
Angstrom, then `atomic_number x y z` per atom. It writes
`from_extpot_to_eon`: one total energy in eV, followed by one force triple
in eV/Angstrom per atom in the same order. The script uses only Python's
standard library. eOn supplies the exchange directory as its working
directory; no absolute path to the source checkout is embedded in it.
