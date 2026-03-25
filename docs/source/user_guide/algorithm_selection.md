---
myst:
  html_meta:
    "description": "Guide to selecting the appropriate algorithm in eOn for your simulation task."
    "keywords": "eOn algorithm selection, NEB vs dimer, optimizer comparison, method guide"
---

# Choosing an algorithm

This page helps you select the right eOn algorithm for your task.

## Finding transition states

| Method | Use when | Configuration |
|---|---|---|
| **NEB** | You know both reactant and product | <project:neb.md> |
| **Dimer** | You have one minimum and want to explore | <project:dimer.md> |
| **Process Search** | Automated saddle point exploration from one state | <project:process_search.md> |
| **Saddle Search** | Single saddle point search with displacement | <project:saddle_search.md> |

### NEB vs Dimer

Use **NEB** when you have both endpoint structures and want the minimum energy
path. NEB optimizes the entire path simultaneously and identifies the
transition state as the highest-energy image.

Use the **Dimer** method when you have only one minimum and want to find
nearby saddle points. The dimer walks uphill along the lowest curvature
direction without needing the product state.

**Process Search** combines both: it uses the dimer to find saddle points,
then minimizes from the saddle to identify the product.

## Accelerated dynamics

| Method | Acceleration type | Use when |
|---|---|---|
| **Parallel Replica** | Spatial parallelism | You have many replicas available | <project:parallel_replica.md> |
| **TAD** | Temperature extrapolation | Barriers are known to follow Arrhenius | <project:saddle_search.md> |
| **Hyperdynamics** | Bias potential | Barriers are localized (bond-boost) | <project:hyperdynamics.md> |
| **AKMC** | Kinetic Monte Carlo | Long timescale evolution with rare events | <project:akmc.md> |

## Optimization

| Optimizer | Convergence | Robustness | Memory |
|---|---|---|---|
| **LBFGS** | Fast near minima | Can fail far from min | O(memory * N) |
| **FIRE** | Moderate | Very robust | O(N) |
| **CG** | Moderate | Good | O(N) |
| **QuickMin** | Slow | Very robust | O(N) |
| **SD** | Very slow | Guaranteed descent | O(N) |

For most tasks, start with **LBFGS**. If it fails to converge, try **FIRE**
or the refinement feature (start with FIRE, switch to LBFGS at threshold).

See <project:optimizer.md> for configuration details.

## Global optimization

| Method | Use when |
|---|---|
| **Basin Hopping** | Finding global minimum of a cluster or surface | <project:basin_hopping.md> |
| **Monte Carlo** | Sampling configurations at finite temperature |
| **Replica Exchange** | Overcoming barriers via temperature exchange |
