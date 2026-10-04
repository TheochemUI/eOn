# Relax engine library

`libeon_relax_engine` lets a host that owns a surface, such as a
Gaussian-process model, run eOn's nudged elastic band or saddle search
against that surface. The host supplies energies, forces and variances
for a batch of images. eOn owns the band, the optimizer and the
convergence logic.

The C header is `include/eon/relax/eon_relax_engine.h`. It is the
reverse of `engine_c_abi.h` used for rgpot engines: there eOn calls the
engine, here the host calls eOn.

| Symbol | Purpose |
|---|---|
| `eon_relax_create` | Build an engine from a Cap'n Proto `RelaxEngineParams`. A NULL config selects an NEB with the `Parameters` defaults. |
| `eon_relax_run` | Run the band or the saddle search to its stopping rule. |
| `eon_relax_step`, `eon_relax_reset` | Advance one optimizer step, or restart the stepper from new endpoints. |
| `eon_relax_destroy` | Release the engine. |
| `eon_relax_version_hash_str` | Identity of the build: the version plus the full git hash. |

`RelaxEngineParams` lives in `schema/eon_relax_engine.capnp`. Its `kind`
field is `neb` or `saddle`, and `NebParams` and `SaddleParams` carry
the image count, iteration cap and force tolerance. `surfaceEpoch`
tells the engine which generation of the host surface the call belongs
to.

Return values are three-valued. Zero means the call completed and the
outcome struct holds the result. A positive value is recoverable. A
negative value is a named failure, and only `EON_RELAX_SURFACE_FATAL`
leaves the engine unusable.

## Surface epochs

A refit surface returns different energies at identical positions.
`Matter::setSurfaceEpoch` and `Potential::surfaceEpoch` key the energy
and variance caches on the surface generation, so a refit invalidates
them without a change in positions.

## Build and test

The library builds with the default configuration and installs next to
`libeonclib`. The unit test `test_relax_engine` links the library and
also opens it with `dlopen` through the `EON_RELAX_ENGINE` environment
variable that `meson test` sets:

```bash
meson test -C bbdir test_relax_engine --print-errorlogs
```
