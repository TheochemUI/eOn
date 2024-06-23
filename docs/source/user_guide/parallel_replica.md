# Parallel Replica

Parallel Replica dynamics (PRD) is the simplest and the accurate way to
accelerate a molecular dynamics simulation as discussed by
{cite:t}`prd-voterParallelReplicaMethod1998` and more recently reviewed by
{cite:t}`prd-perezParallelReplicaDynamics2015`. The only assumption made in this
method is that the reactions satisfy first order kinetics.

```{math}
\mathrm{Pr}(t) = k \exp (-k t)
```

PRD boosts the simulation linearly with the number of replicas and can be easily
combined with other methods for extending the MD time scale, e.g.  the
hyperdynamics method, giving a multiplicative effect in the time scales that can
be achieved.

In the PRD approach, {math}`N` replicas of the system are made at first and then the
momentum in each replica is randomized and dephasing stage is employed to
decorrelate their motions. The simulation clock starts after this dephasing
stage and stops when the first transition is detected in any of replicas.
Because those {math}`N` trajectories are independent, they can explore the phase space {math}`N`
times faster than using a single trajectory. The overall simulation clock is
advanced by the sum of all the simulation times in replicas.

In order to work with distributed computing, we have modified the traditional
scheme for running PRD. The replica generating and dephasing stage is exactly
the same.  However, we make all replicas run the same number of MD steps to
avoid biasing the successful transition trajectories. In other words, results
will only be reported back when the clients finish their full trajectories. The
server increments the simulation time {math}`t` until the first transition
occurs.

In order to run Parallel Replica jobs:
- Set {any}`job <eon.schema.MainConfig.job>` to *parallel_replica* in the
**[Main]** section.
- For regular MD the {any}`time step <eon.schema.DynamicsConfig.time_step>` and
{any}`length <eon.schema.DynamicsConfig.time>` of the trajactory and parameters
of thermostat can be set in the **[Dynamics]** section.
- The {any}`temperature <eon.schema.MainConfig.temperature>` for the dynamics
  run is set in the **[Main]** section.

## Configuration

```{code-block} ini
[Parallel Replica]
```

```{versionchanged} 2.1_TBA
In TOML, this will be `[Parallel_Replica]`
```


```{eval-rst}
.. autopydantic_model:: eon.schema.ParallelReplicaConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: PRD_
keyprefix: prd-
---
```
