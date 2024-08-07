# Matter design

The design of Matter is meant to be orthogonal to the reading and writing of
different outputs. This is conceptually close to the design of `Potential`.

We need to set some boundaries about what Matter is and what it is not meant to
be. In practice, it is desirable to be able to construct empty matter objects,
especially because having to correctly initialize `Matter` before use can be
restrictive, i.e. lazy loading of data is desirable.


```{versionchanged} 3.x
- Matter objects do not hold headers, these are stored in ConFileParser
```

```{versionadded} 3.x
```


- `Matter` is now populated by `ConFileParser`, which handles validation and works for both `con` and `convel` files .
