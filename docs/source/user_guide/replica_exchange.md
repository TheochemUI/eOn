---
myst:
  html_meta:
    "description": "Replica exchange molecular dynamics in the eOn client."
    "keywords": "eOn replica exchange, temperature ladder, Metropolis swap"
---

# Replica exchange

`job = replica_exchange` reads `pos.con`. It writes `results.dat` and the
coldest replica as `pos_out.con`.

`[Replica Exchange] replicas` defaults to 10. `temperature_distribution`
is `exponential` (the default) or `linear`. Any other value stops the
job. `temperature_low` defaults to `[Main] temperature`. A low temperature
that is not positive is replaced by `[Main] temperature`, or by 300 K
when that temperature is not positive either. When `temperature_high`
is not above `temperature_low`, the job sets it to 1.5 times
`temperature_low`.

`sampling_time` (default 1000) and `exchange_period` (default 100) are in
femtoseconds. The integrator step is `[Dynamics] time_step`. When
`thermostat` is omitted, the client uses `andersen`.

At each exchange period the job proposes one swap per replica, between
neighboring temperatures. After the file is read, `exchange_trials` is
set equal to `replicas`, so the key does not change that count. A
proposed swap of energies `E_high` and `E_low` is accepted with
probability

```{math}
\min\left(1, \exp\left[(E_{high}-E_{low})\left(\frac{1}{k_B T_{high}}-\frac{1}{k_B T_{low}}\right)\right]\right).
```

```{code-block} ini
[Main]
job = replica_exchange
temperature = 300

[Replica Exchange]
replicas = 8
temperature_distribution = exponential
temperature_low = 300
temperature_high = 600
sampling_time = 1000
exchange_period = 100

[Dynamics]
time_step = 1.0
thermostat = langevin
```
