# In-process jobs (no file IPC)

The AKMC server used to launch `eonclient` as a subprocess that read
`pos.con` / `config.ini` and wrote `results.dat`. pyeonclient can keep
the Matter in memory and call the same job objects.

## How a job is entered

| Path | Entry | I/O |
|---|---|---|
| Classic `eonclient` | `Job::run()` | CON + `results.dat` |
| pyeonclient / LocalInProcess | `runFromMatter(Matter)` | No CON read |

`runFromMatter` is implemented for ReplicaDynamics, ParallelReplica,
ReplicaExchange, GPSurrogate, SaddleSearch, and ProcessSearch. The
subprocess path still exists for servers that have not switched.

## Potential layout flags

`Potential::layoutFlags()` reports how a pot is executed:

- `InProcess` — force() is a function call in this address space
- `NeedsWorkingDirectory` — the pot reads files from cwd (`in.lammps`, VASP)
- `Subprocess` — the pot is another process (ExtPot, some LAMMPS workers)

In-process jobs should not wrap a `NeedsWorkingDirectory` pot without
setting cwd first.

## What still writes files

`saveData` / `results.dat` still exist so the classic server can parse
them. Filling a JobResult envelope from memory is `eOn-4mrf`.
