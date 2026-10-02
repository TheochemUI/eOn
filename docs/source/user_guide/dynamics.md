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

Six thermostat options are available:

| Thermostat | Key | Description |
|---|---|---|
| **Andersen** | `andersen` | Stochastic velocity reassignment with collision probability per step |
| **Nose-Hoover** | `nose_hoover` | Deterministic extended-system thermostat (chains of length 2) |
| **Langevin** | `langevin` | Stochastic friction + random force, good for non-equilibrium |
| **None** | `none` | NVE ensemble (constant energy, no temperature control) |
| **PILE** | `pile` | Ring-polymer Langevin equation on the normal modes |
| **PIGLET** | `piglet` | Normal-mode GLE on the internal modes, Langevin on the centroid |

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

Langevin and Nose-Hoover treat each Cartesian component on its own. A
fixed component is stored with zero velocity, so the Nose-Hoover step
leaves it in place, and the chain counts only the free components.
Langevin draws no random force on a fixed component and writes no step
there.

## Path-integral sampling

`pile` and `piglet` integrate a ring polymer. The velocity Verlet step
above stays the classical update and does not move the beads. Each step
evaluates every bead in one force batch, so calculator groups carry the
beads together. A fixed component is left out of that update: its force
and its momentum stay zero, and the bead step skips it.

`path_beads` is the bead count. `path_springs = trotter` uses the
primitive ring-polymer frequencies. `path_springs = eco` uses economised
springs fitted to harmonic radii of gyration up to `path_eco_omega_max`.
Economised springs are refused with `piglet`. The instanton refuses them
as well, because that discretisation assumes Trotter springs.

`path_eco_omega_max` is an angular frequency in inverse internal time
units. One internal time unit is 10.18 fs, so 1.0 is 9.82e13 rad/s and
hbar omega = 0.0647 eV.

`piglet` reads one drift matrix and one covariance per internal mode from
`path_gle_file`. The file holds the mode count and the matrix size, then
for each mode the drift matrix and the covariance, row by row. The drift
matrix is in inverse internal time units. The covariance is in kelvin:
kB times it is the covariance of the mass-scaled extended momenta in eV,
so a canonical GLE at the ring temperature has N T on the diagonal, with
N = `path_beads` and T the `[Main]` temperature. A file with one matrix
pair fewer than the bead count covers the internal modes; one with a pair
per bead skips the first, centroid, pair. The centroid keeps a separate
Langevin thermostat.
`path_pile_tau` is that damping time in femtoseconds, and
`path_pile_scale` multiplies the critical damping of the internal modes
when the thermostat is `pile`. `path_seed` seeds the path-integral
random numbers.

The sampler can hold the centroid on a hyperplane and average the
Cartesian force along the normal. That average is the mean force for
thermodynamic integration along a band. Its negative is the derivative
of the potential of mean force.

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
- <project:optimizer.md> for structure optimization
