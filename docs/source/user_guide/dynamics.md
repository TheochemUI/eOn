---
myst:
  html_meta:
    "description": "Configuration for running classical molecular dynamics (MD) simulations in eOn based on Newton's equations of motion."
    "keywords": "eOn molecular dynamics, MD, simulation, thermostat, dynamics"
---

# Dynamics

Molecular dynamics based on Newton's classical equations of motion, integrated
with the velocity Verlet algorithm.

```{note}
For production MD simulations, consider using integrations with LAMMPS or ASE
for more efficient dynamics. The eOn dynamics engine is primarily used as a
building block for accelerated methods (Parallel Replica, TAD, Hyperdynamics).
```

To run a standalone dynamics simulation, set **job** to *dynamics* in the
**[Main]** section:

```{code-block} ini
[Main]
job = dynamics
temperature = 300

[Dynamics]
time_step = 1.0
time = 1000.0
thermostat = andersen
```

Or via [rgpycrumbs](https://rgpycrumbs.rgoswami.me):

```python
from rgpycrumbs.eon.helpers import write_eon_config

config = {
    "Main": {"job": "dynamics", "temperature": 300},
    "Dynamics": {"time_step": 1.0, "time": 1000.0, "thermostat": "andersen"},
}
write_eon_config(config, Path("config.ini"))
```

## Thermostats

Four thermostat options are available:

| Thermostat | Key | Description |
|---|---|---|
| **Andersen** | `andersen` | Stochastic velocity reassignment with collision probability per step |
| **Nose-Hoover** | `nose_hoover` | Deterministic extended-system thermostat (chains of length 2) |
| **Langevin** | `langevin` | Stochastic friction + random force, good for non-equilibrium |
| **None** | `none` | NVE ensemble (constant energy, no temperature control) |

### Andersen thermostat

Controls temperature via random velocity reassignment. The collision period
determines how frequently atoms are thermalized:

```{code-block} ini
[Dynamics]
thermostat = andersen
andersen_alpha = 1.0
andersen_collision_period = 100.0
```

### Langevin thermostat

Applies friction and random forces. The friction coefficient controls the
coupling strength to the heat bath:

```{code-block} ini
[Dynamics]
thermostat = langevin
langevin_friction = 0.01
```

## Time parameters

All times are specified in femtoseconds (fs). The internal time unit conversion
is handled automatically.

- `time_step`: integration timestep (default: 1.0 fs)
- `time`: total simulation time (default: 1000.0 fs)

The number of steps is computed as `floor(time / time_step)`.

## Configuration

```{code-block} ini
[Dynamics]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.DynamicsConfig
```

## See also

- <project:parallel_replica.md> for accelerated dynamics via replica parallelism
- <project:hyperdynamics.md> for bias-potential acceleration
- <project:optimizer.md> for structure optimization (not dynamics)
