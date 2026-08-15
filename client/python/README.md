# pyeonclient

In-process eOn client for Python. Algorithm objects on **Matter**
(not a workdir wrapper).

| Class / function | Role |
|------------------|------|
| `Matter` | Geometry + PEF (ASE `Atoms` analogue) |
| `ImprovedDimer` / `Dimer` / `Lanczos` / `Davidson` | Min-mode |
| `MinModeSaddleSearch` | Single-ended saddle search |
| `NudgedElasticBand` | NEB band (`list[Matter]`) |
| `Hessian` / `get_prefactors` | Vibrational analysis / HTST |
| `make_job` + `Job.run` | Full `JobType` factory (workdir jobs) |

Docs: [Python API](https://eondocs.org/user_guide/pyeonclient.html)

## Install

```bash
pip install pyeonclient
pip install 'pyeonclient[ase]'           # Matter ↔ ASE geometry helpers
pip install 'pyeonclient[models]'        # optional Pydantic DimerSpec / NebSpec
# uv: uv pip install 'pyeonclient[models]'
# Core Dimer/NEB work without pydantic; C++ + light Python still enforce
# accelerant="gp" only with method="improved".
```

## Example (dimer → saddle)

```python
import numpy as np
import pyeonclient as pyec

params = pyec.Parameters()
params.potential = pyec.PotType.LJ
pot = pyec.make_potential(params.potential, params)
matter = pyec.Matter(pot, params)
# ... set positions / cell / numbers ...

matter.relax()
mode0 = np.random.default_rng(0).normal(size=matter.positions.shape)
dimer = pyec.ImprovedDimer(matter, params, pot)
dimer.compute(matter, mode0)

ss = pyec.MinModeSaddleSearch(
    matter, dimer.eigenvector, matter.potential_energy, params, pot
)
status = ss.run()
print(pyec.saddle_status_message(status), ss.eigenvalue)
```

## NEB

```python
path = [pyec.from_ase(img, pot, params) for img in images]
neb = pyec.NudgedElasticBand(path, params, pot)
neb.compute()
# Stamped ConFrames for plots (no neb.con write required):
frames = neb.path_frames()  # list[readcon.ConFrame]
path = list(neb.path_images())
```

## Min / saddle movies (in-memory)

```python
out, ok = matter.relax(retain_frames=True)  # default inplace=False returns working copy
movie = out.movie_frames()  # list[readcon.ConFrame]

ss = pyec.MinModeSaddleSearch(matter, mode, e0, params, pot)
ss.run_retain_frames()
climb = ss.climb_frames()
```

## Build

```bash
meson setup build -Dwith_pyeonclient=true
meson compile -C build
pytest tests/test_pyeonclient_eigenmode.py tests/test_pyeonclient_neb.py -v
```
