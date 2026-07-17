# pyeonclient on PyPI

`pyeonclient` is a **separate** PyPI project from `eon-akmc`. Distribution follows
the **torch model**: wheels and sdist on PyPI; optional fat Metatomic builds
against **PyPI torch**, not a conda-forge gate.

| Project | Import | Contents |
|---------|--------|----------|
| `eon-akmc` | `eon` | Server (AKMC) |
| `pyeonclient` | `pyeonclient` | Client core: Matter, Parameters, Potential, Jobs, NEB, RGPOT |

## Wheel variants

| Variant | Meson | Torch in wheel NEEDED? | Use |
|---------|-------|------------------------|-----|
| **base** (default CI) | `with_rgpot=true`, `with_metatomic=false` | No | LJ/EMT/… + **RGPOT** (dlopen engines) |
| **metatomic** | `with_metatomic=true`, `with_rgpot=true` | Yes (links libtorch stack) | Fat Metatomic pot + engine |

```bash
pip install pyeonclient                 # base
pip install 'pyeonclient[metatomic]'    # runtime torch stack (extras only)
# fat linked extension:
PYEONCLIENT_VARIANT=metatomic ./scripts/pyeonclient_build_wheel.sh
```

Probes::

    import pyeonclient as pyec
    pyec.built_with_rgpot()
    pyec.built_with_metatomic()

## Matrix (GIL + free-threaded)

| Artifact | Tag |
|----------|-----|
| Stable ABI base | `cp312-abi3-manylinux_*` |
| Free-threaded 3.13t | `cp313t-manylinux_*` |
| Free-threaded 3.14t | `cp314t-manylinux_*` |
| Metatomic (optional CI) | same abi3 tag, built with `pip install torch` in manylinux |

## Files

| File | Role |
|------|------|
| `pyproject-pyeonclient.toml` | Package identity, extras, default meson args |
| `scripts/pyeonclient_prepare_pyproject.sh` | Swap root pyproject for builds |
| `scripts/pyeonclient_build_wheel.sh` | Local base/metatomic wheel build |
| `.github/workflows/pyeonclient-wheels.yml` | CI wheels + sdist + publish |

## Release

1. Bump `version` in `pyproject-pyeonclient.toml` and
   `client/python/bind/module.cpp` `__version__` (lockstep). Pin
   `pyeonclient[models]` → `eon-schema` floor if schema public API changed.
2. Tag: `git tag -s pyeonclient-v0.3.2 -m "pyeonclient 0.3.2"`
3. Push tag → `pyeonclient-wheels.yml` builds and publishes (CI publish job).

## Cookbook / RGPOT metatomic

- **Fat**: `built_with_metatomic()` True → `potential = Metatomic`.
- **Thin host**: base wheel + `potential = RGPOT`, `backend = metatomic`,
  `engine_path` / `RGPOT_METATOMIC_ENGINE` → `libmetatomic_engine.so`
  (engine may be produced by a metatomic-variant build).

Neither path requires conda-forge.
