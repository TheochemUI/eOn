---
myst:
  html_meta:
    "description": "Overview of the EON v2.8.0 release, highlighting major features like RO-NEB-CI, metatomic potential support, and a modernized build system."
    "keywords": "EON v2.8.0, release, RO-NEB-CI, metatomic, GPR-Dimer, Meson build"
---

## [v2.8.0] - 2025-09-04

First release in a decade.

This is a landmark release for EON, representing over five years of dedicated
development. The entire framework has been modernized to be more powerful,
easier to use, and more accessible to researchers across all major platforms.
This release introduces cutting-edge simulation methods, a vast expansion of
supported computational chemistry codes, and a complete overhaul of the build
and installation process.

### Highlights

* **Advanced Transition State Finding with RO-NEB-CI**: A powerful new tool for
  accurately locating complex transition states.
* **Machine-Learned Potentials with `metatomic`**: Full support for `metatomic`
  potentials via the [metatensor
    ecosystem](https://docs.metatensor.org/latest/index.html).
* **Expanded Potential Interfaces**: Drastically increased interoperability with
  direct interfaces for **NWChem**, **ORCA**, **AMS**, **XTB**, and **ASE**.
* **Input Validation and Schema**: A comprehensive Pydantic schema for all
  configuration files provides automatic input validation and clearer error
  messages.
* **Modernized Build System**: A complete overhaul using Meson, providing a
  faster, more reliable, and truly cross-platform installation experience.

This release is also supported by several new publications demonstrating the
implementation and application of these new methods. For a complete record of
all changes and accompanying research, please explore the detailed pages for
this version:

```{toctree}
:maxdepth: 2
:caption: Contents

release-notes
publications
```
