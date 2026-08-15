---
myst:
  html_meta:
    "description": "Run eonclient jobs under AiiDA with the aiida-eon plugin."
    "keywords": "eOn, AiiDA, aiida-eon, CalcJob, workchain, Slurm"
---

<!-- vale proselint.Uncomparables = NO -->
# AiiDA (`aiida-eon`)

[aiida-eon](https://pypi.org/project/aiida-eon/) is the AiiDA 2 plugin for
the `eonclient` binary. One calculation writes `config.ini` and the job's
`.con` files, runs `eonclient`, and parses `results.dat`.

The Python adaptive kinetic Monte Carlo (AKMC) server
(`python -m eon.server`) is a different process.
`WorkflowFactory("eon.akmc")` submits independent `process_search`
clients; superbasin kinetic Monte Carlo still belongs on the server.
In-process Matter work is [pyeonclient](pyeonclient.md).

```{seealso}
Plugin repository: [HaoZeke/aiida-eon](https://github.com/HaoZeke/aiida-eon).
`config.ini` authorship is {func}`~rgpycrumbs.eon.helpers.write_eon_config`
({doc}`/tutorials/dict_config`).
```

## Install

```{code-block} bash
pip install aiida-eon
verdi plugin list aiida.calculations eon
```

Register a prebuilt `eonclient`:

```{code-block} bash
verdi code create core.code.installed \
  --label eonclient \
  --computer localhost \
  --filepath-executable "$(command -v eonclient)" \
  --default-calc-job-plugin eon
```

On a Slurm login node the computer is `core.local` plus `core.slurm`.
The repository includes Elja YAML under `examples/elja/`.

## A minimization

<!-- vale off -->
```{code-block} python
from aiida import orm
from aiida.engine import run
from aiida.plugins import CalculationFactory, DataFactory

EonCalculation = CalculationFactory("eon")
ConData = DataFactory("eon.con")
result = run(
    EonCalculation,
    code=orm.load_code("eonclient@localhost"),
    parameters=orm.Dict(
        dict={
            "Main": {"job": "minimization"},
            "Potential": {"potential": "lj"},
            "Optimizer": {"opt_method": "lbfgs", "converged_force": 0.01},
        }
    ),
    structure=ConData.from_path("pos.con"),
)
print(result["results"].get_dict()["potential_energy"])
```
<!-- vale on -->

The `job` key selects which input files the client expects.
Aliases such as `neb` become `nudged_elastic_band` before `config.ini`
is written. A nudged elastic band (NEB) needs `reactant` and
`product`. Extra potfiles (LAMMPS, Vienna Ab initio Simulation Package)
go in the `potfiles` `FolderData` port and are copied into the client
working directory.

## Workchains

```{code-block} python
from aiida.plugins import WorkflowFactory

Minimize = WorkflowFactory("eon.minimize")
Neb = WorkflowFactory("eon.neb")
Saddle = WorkflowFactory("eon.saddle")
Process = WorkflowFactory("eon.process_search")
Prefactor = WorkflowFactory("eon.prefactor")
Akmc = WorkflowFactory("eon.akmc")
```

Each workchain pins the `job` key and exposes the CalcJob ports under
`calc`. Serve mode (`eonclient --serve-*`) and `forceBatch` are not
jobs.

## Related surfaces

| Token | Where it lives |
|---|---|
| `akmc`, `escape_rate` | Python server, or `eon.akmc` for search campaigns |
| `eonclient --serve-*` | Long-lived potential remote procedure call |
| `forceBatch` | Engine path inside NEB and dimer |
| Matter / `relax` / `NEB.compute` | [pyeonclient](pyeonclient.md) |
