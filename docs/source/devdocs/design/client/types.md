# Type system

`eON` was written to use `Eigen` as the matrix multiplication library.
*However*, we define `EIGEN_DEFAULT_TO_ROW_MAJOR` which makes the Eigen
constructs of `eON` not directly interoperable with those from other libraries.

```{versionchanged} 2.x
Generic type interfaces were setup, to eventually transition away from the preprocessor macro.
```

## Header conventions

```{versionchanged} 2.x
```

- `include-what-you-use`
- `pragma once` is preferred over setting `ifdef`s
