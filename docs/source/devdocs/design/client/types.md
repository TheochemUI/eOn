---
myst:
  html_meta:
    "description": "Information for developers on the type system used in the eOn C++ client, particularly the use of the Eigen library."
    "keywords": "eOn type system, C++, Eigen library, row major, header conventions"
---

# Type system

`eOn` was written to use `Eigen` as the matrix multiplication library.
Row-major storage is specified via the `eOnStorageOrder` constant in
`client/Eigen.h`, which the `AtomMatrix` and `RotationMatrix` typedefs use.
This keeps the layout explicit and in one place, without the old
`EIGEN_DEFAULT_TO_ROW_MAJOR` preprocessor macro that silently changed every
Eigen type.

```{versionchanged} 2.x
Replaced the `EIGEN_DEFAULT_TO_ROW_MAJOR` preprocessor macro with an explicit
`eOnStorageOrder` constant used by the project typedefs.
```

## Header conventions

```{versionchanged} 2.x
```

- `include-what-you-use`
- `pragma once` is preferred over setting `ifdef`s
