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
