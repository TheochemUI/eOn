# Matter design

The design of Matter is meant to be orthogonal to the reading and writing of
different outputs. This is conceptually close to the design of `Potential`.

We need to set some boundaries about what Matter is and what it is not meant to
be. In practice, it is desirable to be able to construct empty matter objects,
especially because having to correctly initialize `Matter` before use can be
restrictive, i.e. lazy loading of data is desirable.

## Populating

```{versionchanged} 3.x
- Matter objects do not hold headers, these are stored in ConFileParser
```


```{versionadded} 3.x
```

- `Matter` is now populated by `ConFileParser`, which handles validation and works for both `con` and `convel` files .

## Caching

```{versionadded} 3.x
```

Earlier versions of `eOn` used the concept of a tainted bit to provide a way to
short-circuit recalculations at the same state. Essentially, the
`recomputePotential` boolean was set for certain operations and controlled when
/ how the potential energy and forces were computed. This was essentially an LRU
cache of size one, which notably prevented const-correctness of methods, and
prevented re-use of previous calculations.

```{todo}
```
The new design requires `Matter` to be constructed with a `cachelot` instance, with the parameters set at runtime but offers the ability to re-use older calculations and eventually:

- Run in an async thread (TODO(rg))
- Be used for structure comparison (TODO(rg))
- Allow on-the-fly training of potentials for positions close to existing ones (TODO(rg))
