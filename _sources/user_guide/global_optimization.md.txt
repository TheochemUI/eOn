---
myst:
  html_meta:
    "description": "Minima hopping in the eOn client, under job global_optimization."
    "keywords": "eOn minima hopping, global optimization, molecular dynamics escape"
---

# Minima hopping

`job = global_optimization` is minima hopping. It is a different job from
basin hopping. The client logs `Beginning minima hopping` and reads
`pos.con`. It writes `monitoring.dat` and `earr.dat` in the working
directory.

`[Global Optimization] move_method` is `md` (the default) or `random`.
The molecular-dynamics move is a constant-energy trajectory. It stops
after `mdmin` potential-energy minima (default 3) or after 1000 steps.
The random move uses `[Basin Hopping] displacement` (default 0.5, in
the length unit of the positions) and that section's
`displacement_distribution`.

`decision_method` is the non-probabilistic energy window (`npew`, the
default) or `boltzmann`. Any other value stops the client while the
file is read. `npew` accepts a hop
whose energy is below the current energy plus a window. The window
starts at 0.1 eV. An acceptance multiplies it by `1/alpha`, and a
rejection multiplies it by `alpha`. `alpha` defaults to 1.02.
`boltzmann` keeps an uphill hop with probability
`exp(-de/(kB*temperature))` at `[Main] temperature`. A temperature or
`kB` that is not positive rejects that hop.

`beta` (default 1.05) scales the kinetic energy of the next
molecular-dynamics escape. A repeat of the current minimum, or a minimum
already in the list, multiplies it by `beta`. A new minimum multiplies
it by `1/beta`. `steps` (default 10000) is the number of hops. The loop
also stops when the current energy falls below `target_energy`.

```{code-block} ini
[Main]
job = global_optimization
temperature = 300

[Global Optimization]
move_method = md
decision_method = npew
steps = 10000
beta = 1.05
alpha = 1.02
mdmin = 3
```
