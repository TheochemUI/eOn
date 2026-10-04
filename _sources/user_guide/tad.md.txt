---
myst:
  html_meta:
    "description": "The client job tad runs a hot trajectory and a cold-temperature correction."
    "keywords": "eOn temperature accelerated dynamics, Arrhenius, low temperature"
---

# Temperature accelerated dynamics

Temperature accelerated dynamics runs a hot trajectory at `[Main] temperature`.
The client job is `tad`.
The cold temperature is `[TAD] low_temperature`, default 300 K.
`min_prefactor` (default 0.001) and `confidence` (default 0.001) set how
long that hot trajectory continues after a transition. The trajectory
length and the step come from `[Dynamics]`. When `thermostat` is
omitted, the client uses `andersen`. The trajectory tests for a new
state on the parallel-replica state-check interval. The record interval
in that section can refine the crossing.

A barrier `Eb` on the trajectory scales the hot transition time by
`exp(Eb/kB * (1/low_temperature - 1/temperature))`. `tmin` is the
shortest corrected time. `factor` is `ln(1/confidence) / min_prefactor`.
The hot trajectory can stop before `[Dynamics] time`. It stops once the
hot clock passes `factor * (tmin / factor)` raised to
`low_temperature / temperature`.

```{code-block} ini
[Main]
job = tad
temperature = 600

[TAD]
low_temperature = 300
min_prefactor = 0.001
confidence = 0.001

[Dynamics]
time_step = 1.0
time = 1000.0
```
