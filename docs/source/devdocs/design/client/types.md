---
myst:
  html_meta:
    "description": "Information for developers on the type system used in the EON C++ client, particularly the use of the Eigen library."
    "keywords": "EON type system, C++, Eigen library, row major, header conventions"
---

# Type system

`eOn` was written to use `Eigen` as the matrix multiplication library.
*However*, we define `EIGEN_DEFAULT_TO_ROW_MAJOR` which makes the Eigen
constructs of `eOn` not directly interoperable with those from other libraries.

```{versionchanged} 2.x
Generic type interfaces were setup, to eventually transition away from the preprocessor macro.
```

## Header conventions

```{versionchanged} 2.x
```

- `include-what-you-use`
- `pragma once` is preferred over setting `ifdef`s
