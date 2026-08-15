---
myst:
  html_meta:
    "description": "User guide for the Adaptive Kinetic Monte Carlo (aKMC) method in eOn. Learn how to coarse-grain molecular dynamics for rare event systems."
    "keywords": "aKMC, adaptive kinetic monte carlo, eOn, rare event systems, transition state theory"
---

# Adaptive Kinetic Monte Carlo

The adaptive kinetic Monte Carlo (aKMC) method is a method to coarse grain
molecular dynamics for rare event systems as described in
{cite:t}`ak-henkelmanLongTimeScale2001,ak-xuAdaptiveKineticMonte2008`.  A rare
event system is one in which the interesting dynamics is governed by short
transitions between stable states. The fast vibrational motion within a stable
state is considered to be in equilibrium and described statisitcally. A
transition between states is treated as first order, since it is a rare
event, and the rate of the transition is calculated from the harmonic
approximation to transition state theory (hTST).

The hTST approximation of a transition rate is calculated from the energy
difference between the saddle point along the minimum energy path for the
transition and the initial minimum. The vibrational modes at these points are
also used to calculate the prefactor. An hTST rate is of the standard Arrhenius
form {math}`R = v \exp (-\Delta E/kT)` where {math}`v` is the product of all positive
modes at the minimum divided by those at the saddle, {math}`\Delta E` is the energy
barrier, and {math}`kT` is the thermal enregy.

To propagate the dynamics within aKMC, a list of all possible rates
leading away from the current stable state to any other state is required.
Formally, there are an enourmously large number of such transitions (also called
processes) available in a typical atomic system, but in fact only the
transitions with the fastest rates with the highest probability of happening are
required. The search for processes is then limited to those with rates on the
same order as the fastest processes found.

The search for possible processes is the primary task of the aKMC simulations.
Each client does a minimum mode (minmode) following search from the miminum of
the current state and searches for a saddle point which connects from the
minimum in the current state to an adjacent state
{cite:p}`ak-henkelmanDimerMethodFinding1999`. A saddle point is connected to a
state if a minimization initiated along the negative mode at the saddle
converges to the minimum of that state.

Each client is tasked with one or more such searches. It climbs from the minimum
to a saddle, and if successful, it minimizes on either side of the saddle to
determine the connecting states. The prefactor for the transition is also
calculated by finite difference and the hTST rate is calculated. These data are
reported back to the server.

Adaptive (a), or off-lattice KMC {cite:p}`ak-trochetOffLatticeKineticMonte2020`
works in continuous configuration space. The method builds a rate table from
local atomic movements and energy barriers, so atoms can sit off a grid. The
confidence scheme and thermal accessibility settings restrict accepted
transitions to those that meet the documented thresholds.

The server is reponsible for the time evolution of the system by the KMC
algorithm. Each process leading to a new state is tabulated in a rate table and
one transition is selected stochastically with a probability proportional to its
rate. The transition time is drawn from a first-order distribution for the total
rate of escape from the state.

A complete tutorial is [also provided](project:../tutorials/akmc.md).

## Configuration

```{code-block} ini
[akmc]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.AKMCConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: AK_
keyprefix: ak-
---
```
