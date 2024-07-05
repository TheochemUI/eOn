# Memory management

We follow the [GSL resource management
rules](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#S-resource).
Notably:
- Raw pointers are non-owning
- Raw references are non-owning

```{versionchanged} 2.x
Heap allocations are to be much reduced.
```
