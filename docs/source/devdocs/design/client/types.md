---
myst:
  html_meta:
    "description": "Information for developers on the type system used in the eOn C++ client, particularly the use of the Eigen library."
    "keywords": "eOn type system, C++, Eigen library, row major, header conventions"
---

# Type system

`eOn` was written to use `Eigen` as the matrix multiplication library.
All matrices are row-major (`eOnStorageOrder = Eigen::RowMajor`), so that
`.data()` yields `[x0, y0, z0, x1, y1, z1, ...]` as the Fortran potentials
and `VectorXd::Map` round-trips expect.

Instead of the `EIGEN_DEFAULT_TO_ROW_MAJOR` preprocessor macro (which would
make eOn's Eigen types binary-incompatible with other libraries), `client/Eigen.h`
provides explicit row-major type aliases (`MatrixXd`, `Matrix3d`, `Matrix4d`,
`AtomMatrix`, `RotationMatrix`) alongside selective `using` declarations for
vector types (`VectorXd`, `VectorXi`, `Vector3d`) that are unaffected by
storage order.

```{warning}
New multi-dimensional Eigen matrix types in the codebase **must** use
`eOnStorageOrder` or one of the existing row-major aliases.  Adding a bare
`Eigen::MatrixXd` (column-major) will silently corrupt force projections and
data mapping.
```

```{versionchanged} 2.11
Removed the `EIGEN_DEFAULT_TO_ROW_MAJOR` preprocessor macro.  All matrix types
are now explicitly row-major via aliases in `client/Eigen.h`, ensuring binary
compatibility with other Eigen-based libraries.
```

## Header conventions

```{versionchanged} 2.x
```

- `include-what-you-use`
- `pragma once` is preferred over setting `ifdef`s
