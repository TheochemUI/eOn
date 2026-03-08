---
myst:
  html_meta:
    "description": "Detailed release notes for eOn v2.12.0, covering quill migration, namespace eonc, PotRegistry, SafeMath, and FPE hardening."
    "keywords": "eOn release notes, quill, namespace eonc, PotRegistry, SafeMath, FPE"
---

# Release notes

## [v2.12.0] - 2026-03-08

### Removed

#### _potcalls.log text logger

The per-instance `FileScoped` text logger that emitted `QUILL_LOG_TRACE_L3` on
every `get_ef()` call has been removed. The dead static counters
`Potential::fcalls`, `fcallsTotal`, `wu_fcallsTotal`, and `totalUserTime` are
also removed. Output is now `_potcalls.json` with structured per-instance
records (see PotRegistry below).

### Added

#### PotRegistry

A thread-safe, enum-indexed singleton (`PotRegistry`) tracks per-instance
potential lifecycle: creation time, destruction time, force call count, and a
unique ID. Per-type aggregate counters use relaxed atomics for lock-free
`force_calls`, `created`, and `alive` tracking. Per-instance records are
mutex-protected and appended on destruction.

Force call delta tracking is restored in SafeHyperJob, TADJob,
ParallelReplicaJob, ReplicaExchangeJob, NudgedElasticBandJob, BasinHoppingJob,
HessianJob, PrefactorJob, and MonteCarloJob using
`PotRegistry::get().total_force_calls()`.

Unit tests: 6 test cases, 38 assertions covering lifecycle, force calls, delta
tracking, reset, JSON output, multi-type, and instance records.

#### SafeMath.h

Guarded arithmetic functions (`safe_div`, `safe_recip`, `safe_acos`,
`safe_sqrt`, `safe_atan_ratio`) and an Eigen-aware `safe_normalized` template.
These prevent SIGFPE from division-by-zero and domain errors without changing
results for valid inputs.

#### Quill logging backend

Replaced spdlog/fmt with quill for lock-free asynchronous logging. The new
`EonLogger.h` API provides `eonc::log::Scoped` (RAII console logger) and
`eonc::log::FileScoped` (RAII file logger) helpers, reducing boilerplate from
14 lines to 1.

### Changed

#### namespace eonc

All client classes, enums, and helper namespaces are wrapped under
`namespace eonc`. The `helper_functions` namespace is renamed to `helpers`.
Backward-compatible `using` aliases are provided at global scope. `using
namespace std` is removed from all client headers.

#### Other changes

- Replace cxxopts with argum for CLI parsing (C++20).
  ([#320](https://github.com/TheochemUI/eOn/issues/320))
- Replace per-site `scoped_interpreter` guards with lazy singleton
  `eonc::ensure_interpreter()` in `PyGuard.h`.
  ([#324](https://github.com/TheochemUI/eOn/issues/324))
- Updated metatomic ecosystem: torch 2.10, metatomic-torch 0.1.9+,
  metatensor-torch 0.8.4+, vesin/vesin-torch 0.5.2+, metatrain 2026.2.1+.
- Updated vesin to v0.5.2: per-dimension periodicity (`bool[3]`) and
  VesinDevice struct syntax.
- Switched all logging macros from bare `LOG_*` to `QUILL_LOG_*` with
  `QUILL_DISABLE_NON_PREFIXED_MACROS` to prevent macro collisions.
- Optimized quill backend: reduced sleep duration (100us -> 10us), larger
  transit buffer (256 -> 2048), faster flush interval (200ms -> 100ms).

### Fixed

- Fix ASE_POT compilation errors and rename `-DASE_POT` to `-DWITH_ASE_POT`.
  ([#321](https://github.com/TheochemUI/eOn/issues/321))
- Fix `[IDimerRot]` column misalignment.
  ([#322](https://github.com/TheochemUI/eOn/issues/322))
- Suppress FPE trapping during libtorch operations in MetatomicPotential.
  ([#323](https://github.com/TheochemUI/eOn/issues/323))
- Change `uncertainty_threshold` default from `0.1` to `-1` (disabled).
  ([#325](https://github.com/TheochemUI/eOn/issues/325))
- Fixed quill migration test failures: logger initialization, ConfigParser
  defaults, MSVC `-DNOMINMAX`.
  ([#327](https://github.com/TheochemUI/eOn/issues/327))
- Guard unprotected floating-point divisions and `.normalize()` calls across
  12 source files using `eonc::safemath` utilities.
- Make POSIX FPE signal handler async-signal-safe.
- Use `-isystem` for pip-installed metatomic/vesin includes to suppress
  third-party warnings.
