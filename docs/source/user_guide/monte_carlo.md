---
myst:
  html_meta:
    "description": "Canonical Monte Carlo sampling in the eOn client."
    "keywords": "eOn Monte Carlo, Metropolis, step size, temperature"
---

# Monte Carlo

`job = monte_carlo` reads `pos.con`. When `[Main] checkpoint` is set and
`pos_cp.con` exists, the client reads that file. It writes `out.con`,
`results.dat`, and `movie.con`.

`[Monte Carlo] step_size` (default 0.005) is the standard deviation of a
Gaussian kick on each Cartesian coordinate, in the length unit of the
positions. Fully fixed atoms get no kick. `[Monte Carlo] steps`
(default 1000) is the number
of trials. Both values must be positive.

The temperature is `[Main] temperature`, in kelvin, and must be positive.
A downhill trial is kept. An uphill trial is kept with probability
`exp(-de/(kB*temperature))`. `kB` is in eV/K.

```{code-block} ini
[Main]
job = monte_carlo
temperature = 300

[Monte Carlo]
step_size = 0.005
steps = 1000
```
