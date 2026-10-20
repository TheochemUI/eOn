---
myst:
  html_meta:
    "description": "Guide to coarse-graining methods in eOn, such as MCAMC and AS-KMC, for accelerating simulations with vastly different timescales."
    "keywords": "eOn coarse graining, MCAMC, AS-KMC, superbasins, accelerated simulation"
---

# Coarse graining

In aKMC simulations where there are vastly different rates, the simulation can
get stuck in a group of states connected by relatively fast rates. Exploring
slower transitions can take a prohibitively large number of KMC steps.
`eOn` implements two methods that skip those fast recrossings.

```{note}
AS-KMC and MCAMC cannot be used simultaneously.
```

## Monte Carlo with Absorbing Markov Chains (MCAMC)

The first method, projective dynamics, described in
{cite:t}`cg-novotnyTutorialAdvancedDynamic2001`, groups states that are joined
by fast rates into "superbasins". Information about transitions between states
in a superbasin is lost, but the rates for transitions across a superbasin are
correct.

## Accelerated Superbasin Kinetic Monte Carlo (AS-KMC)

The second method, accelerated superbasin kinetic Monte Carlo (AS-KMC) of
{cite:t}`cg-voterHyperdynamicsAcceleratedMolecular1997`, artificially raises low
barriers. The dynamics between states connected by fast rates are simulated, but
an error is introduced in the dynamics direction and time.


The basic process of AS-KMC involves gradually raising process barriers found to
be inside of a superbasin such that exiting from the basin gradually becomes
more likely. The method is designed to raise all the barriers in the superbasin
simultaneously. Once a particular barrier has been crossed a certain number of
times, {math}`N_f` (more on determining {math}`N_f` shortly), a check is
performed to determine whether or not the current state is part of a superbasin.
This is called the Superbasin Criterion.

In the Superbasin Criterion, a search is performed, originating at the current
state and proceeding outward through all low-barrier processes to adjacent
states, and then through all low-barrier processes from each of these states,
etc. For each low-barrier process found, if the process has been followed fewer
than {math}`N_f` times, the Superbasin Criterion fails and no barriers are
raised.

Thus, in the outward-expanding search from the originating state, the search
continues until either a low-barrier process has been seen fewer than
{math}`N_f` times (and the Criterion fails) or until all connected low-barrier
processes have been found and have been crossed at least {math}`N_f` times (the
edges of the superbasin are then defined and the Criterion passes). If the
Superbasin Criterion passes, all the low-barrier processes (each of which as
been crossed {math}`N_f` times) are raised.

Several parameters dictate the functioning of the AS-KMC method. These
parameters dictate how much the barriers are raised each time the Superbasin
Criterion
passes({any}`eon.schema.CoarseGrainingConfig.askmc_barrier_raise_param`), what
defines a "low-barrier" for use in the Superbasin
Criterion({any}`eon.schema.CoarseGrainingConfig.askmc_high_barrier_def`), and
the approximate amount of error the user might expect in eventual superbasin
exit direction and time compared to normal KMC simulation
({any}`eon.schema.CoarseGrainingConfig.askmc_confidence`).

## amsel discover_decide

`[amsel] discover_decide = true` runs on the current state's process table.
The state list can hold one state. `use_mcamc` stays off.
`use_mcamc` is not required. The repeat-count
confidence scheme does not hold this step. A run that leaves
`discover_decide` off still waits for `confidence` under
`confidence_scheme`.

A barrier strictly below `e_min_init` counts as an in-basin edge. A barrier
at or above `e_min_init` counts as an exit. On Si6N8 isomer 1 the back
barrier is 0.20 eV and the flip-out barrier is 0.30 eV. With `e_min_init`
at 0.25 eV the 0.20 eV edge stays in the basin and the 0.30 eV edge leaves.
A table that has only the 0.30 eV saddle, and no faster edge, still leaves.
The transient set is that state. The product is absorbing. A product column
of -1 is that same exit: the label sent to amsel is a 32-bit id at or above
2147483648, and the hop creates the product state from the process id.

The basin comes from `amsel.discover_decide_status` when a faster edge is
present. The exit time and the exit channel come from the mean-rate method
(MRM) or from first-passage-time analysis (FPTA). `debug_use_mean_time`
selects MRM. The mean exit time is the MRM value `tau_total`. Otherwise,
FPTA draws one first-passage time. Both kernels live in the `amsel`
package. With that package absent, the log line reads
`amsel discover_decide status=unavailable` and no step is taken until the
confidence threshold is met. Once that threshold is met, the step is
ordinary kinetic Monte Carlo.

```{code-block} ini
[amsel]
discover_decide = true
e_min_init = 0.25
e_min_step = 0.05
e_min_floor = 0.05
cv_threshold = 10.0
on_error = fallback_single
```

One status list:

- `unavailable`, `available` off: `amsel` did not import, or `on_error` is `unavailable_mcamc` after a failed call
- `fallback_single`, `available` on: the call failed and `on_error` is `fallback_single`; the step stays ordinary kinetic Monte Carlo
- `accepted` or `retightened`: one basin, and the exit comes from MRM or FPTA
- `split_required`: the exit is the primary basin that contains the entry state
- `rejected_no_metastable_basin`: no basin, and the step stays ordinary kinetic Monte Carlo

`unavailable` together with `available` on sits outside this list.

```{eval-rst}
.. autopydantic_model:: eon.schema.AmselConfig
```

With `use_mcamc` on as well, an accepted basin still leaves through MRM or
FPTA. A rejected basin falls through to the MCAMC superbasin step.
`split_required` on that MCAMC object keeps the entry state and writes the
reduced member list.

## Configuration

```{code-block} ini
[Coarse Graining]
```


```{eval-rst}
.. autopydantic_model:: eon.schema.CoarseGrainingConfig
```


## References


```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: CG_
keyprefix: cg-
---
```
