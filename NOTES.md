# Bugs found

- `client/Potential.cpp:458` sent every `potential = lammps` configuration
  directly to `LAMMPSPot`, so the existing rgpot integration could never
  select native tabulated EAM. `client/Potential.cpp:172-256` now recognizes
  the safe single-element `eam`, `eam/alloy`, and `eam/fs` subset, chooses
  `EAMTablePot`, and retains direct liblammps as the generic fallback.
- `client/Potential.cpp:583` constructed `RgpotPot` for every RGPOT backend.
  Consequently `backend = lammps` reached a dispatcher that did not implement
  it. The factory now maps `lammps`, `lammpc`, and `lammpsc` to eOn's direct
  in-process LAMMPS library interface.
- `client/potentials/Rgpot/RGPotEngine.cpp:249-333` listed only five engine
  families. It now exposes native EAMTable and every native name built by the
  pinned rgpot profile: LJ, Morse, Sutton-Chen, published Stillinger-Weber,
  Tersoff Si-C, EMT, cell-list EAM, and the enabled Fortran kernels.
- The previous rgpot pin did not contain native tabulated EAM.
  `subprojects/rgpot.wrap:5` now pins the current rgpot develop revision,
  and `client/meson.build:556-587` requires rgpot 3.7 and enables
  `with_eam_table`. No rgpot source change is needed.
- The pinned rgpot's `LammpsPot` is an external-engine conformance fixture for
  a parameterized Lennard-Jones calculation. It is not a general LAMMPS input
  interpreter. Treating it as the Mishin fallback would silently run the wrong
  potential, so eOn's linked liblammps interface is the correct generic
  in-process fallback.
- `client/potentials/LAMMPS/LAMMPSPot.cpp:45-87` enabled the POSIX worker for
  every production instance. `client/potentials/LAMMPS/LAMMPSPot.cpp:418-442`
  now evaluates in-process by default and retains the worker only as an
  explicit compatibility and benchmark mode.
- `client/potentials/LAMMPS/LAMMPSPot.cpp:809-849` previously recognized real
  units only through a private comment marker. A normal `units real` or
  `units metal` line was loaded after `create_box`, where LAMMPS rejects a
  unit change. eOn now reads and applies the supported unit command before
  creating the box and does not replay it.
- `client/potentials/LAMMPS/LAMMPSPot.cpp:418` discarded the variance pointer
  instead of initializing its value. The implementation now writes the
  deterministic-potential variance explicitly.
- `eon/communicator_inprocess.py:452-526` rebuilt Parameters, Potential, and
  Matter for each submitted search. It also preferred invariant configuration
  even when a job supplied its own `config.ini`. It now gives job-local
  configuration precedence and caches one potential and minimized reactant per
  state, invalidating the entry when geometry or configuration changes.
- `eon/config.yaml:1548` and `meson_options.txt:59` left the subprocess
  communicator and disabled Python binding as defaults. The in-process
  communicator and its binding are now the defaults, while explicit
  communicator selection remains available.
- `client/potentials/LAMMPS/LAMMPSPot.cpp:93-177` replayed four `unfix`
  commands, four group deletions, and every fixed-coordinate group on every
  force call. Matter supplies the same all-free mask on each call, so the
  worker spent most of its live AKMC time processing commands that could not
  affect a force. All-free masks are now normalized to empty, unchanged masks
  are ignored, and a changed mask is applied once per LAMMPS instance.
- `eon/communicator_inprocess.py:517-540` ignored `path_pot`. The local
  communicator copied `in.lammps` and its table into each job directory, but
  the in-process path constructed its potential in the state directory. That
  made native EAM detection fail and could fall through to liblammps without
  an input file. Potential support files now have one persistent,
  traversal-checked workspace per communicator.
- `eon/communicator_inprocess.py:57-64,551-589` included the per-search random
  seed in the state-session cache key. Every search therefore invalidated the
  supposedly persistent potential and minimized reactant. Random seeds are
  now excluded from the session identity; geometry and all behavior-affecting
  configuration remain in it. `eon/config.py:233-240` also validated
  `reuse_state_sessions` but never copied its parsed value onto ConfigClass,
  making the comparison switch always use its constructor default. The
  runtime configuration now carries the parsed boolean.
- `eon/explorer.py:140-145` assumed every communicator completed
  asynchronously. An in-process result stayed only in memory, while the AKMC
  invocation exited before its next registration pass and discarded it.
  Synchronous communicators now register their completed batch immediately.
- `client/ProcessSearchJob.cpp:192-252` regenerated a displacement when called
  through the Python binding instead of consuming the server's exact
  `displacement.con` and `direction.dat`. The in-process entry now accepts
  those objects directly, preserves reactant atom indexing and fixed
  coordinates, and follows the same search as the file path.
- `client/ProcessSearchJob.cpp:276-282` did not expose its final status and
  always wrote legacy result files from the in-process entry.
  `eon/communicator_inprocess.py:390-438` also emitted a mapping with keys that
  did not match the typed JobResult schema. The binding now returns status,
  force-call phases, barriers, energies, and prefactors without writing job
  files, and the communicator publishes a complete typed and legacy view.
- `client/python/pyeonclient/__init__.py:99` did not export the already-bound
  `ProcessSearchJob`, so Python silently used the lower-level saddle-only
  fallback. The fallback omitted endpoint minimization, connectivity,
  prefactors, and process status. The complete process-search class is now
  the public in-process path.
