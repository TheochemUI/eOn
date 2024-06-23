# Coarse Graining

In aKMC simulations where there are vastly different rates, the simulation can
get stuck in a group of states connected by relatively fast rates. In order to
explore slower transitions, a prohibitively large number of KMC steps may be
needed. In order to circumvent this problem, `eOn` implements two methods. 

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
defines a “low-barrier” for use in the Superbasin
Criterion({any}`eon.schema.CoarseGrainingConfig.askmc_high_barrier_def`), and
the approximate amount of error the user might expect in eventual superbasin
exit direction and time compared to normal KMC simulation
({any}`eon.schema.CoarseGrainingConfig.askmc_confidence`).

## Configuration

```{code-block} ini
[Coarse Graining]
```

```{versionchanged} 2.1_TBA
In TOML, this will be `[Coarse_Graining]`
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
