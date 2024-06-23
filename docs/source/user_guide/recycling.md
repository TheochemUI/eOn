# Recycling

As reported in {cite:t}`recyc-xuAdaptiveKineticMonte2008`, `eON` implements a
method of saddle point recycling that can significantly reduce the computational
cost of the aKMC algorithm.

Suppose we are in reactant state {math}`R_0`, and we have discovered a series of
saddles and their corresponding products, {math}`S_i` and {math}`P_i`,
respectively. Once we have reached confidence that we have found all
energetically relevant processes, we select one of these processes and move to
the corresponding product state.

For this example, let us assume that we have selected the process with saddle
{math}`S_0` and product {math}`P_0`. If we have found {math}`N` energetically
relevant processes in state {math}`R_0`, we can make suggestions of the saddle
geometries {math}`G_i` for saddles leading out of state {math}`P_0`, i.e.:

```{math}
G_i = P_0 + (S_i - R_0)
```


We use a min-mode following algorithm to converge these suggested saddle points.
To reach confidence again in state {math}`P_0`, we need only perform saddle
searches in the region around the atoms that moved significantly from state
{math}`R_0` to state {math}`P_0`, resulting in a significant reduction in
computational costs. If this region is local, the overall cost does not
increase with the total system size.

## Configuration

```{code-block} ini
[Recycling]
```


```{eval-rst}
.. autopydantic_model:: eon.schema.RecyclingConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: RECYC_
keyprefix: recyc-
---
```
