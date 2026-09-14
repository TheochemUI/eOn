---
myst:
  html_meta:
    "description": "eOn v2.8.0: RO-NEB-CI, metatomic potential support, and a Meson build."
    "keywords": "eOn v2.8.0, release, RO-NEB-CI, metatomic, GPR-Dimer, Meson build"
---

## [v2.8.0] - 2025-09-04

First release in a decade.

Five years of work: new transition-state methods, more potential
interfaces, and a Meson build that runs on Linux, macOS, and Windows.

### Highlights

* **RO-NEB-CI**: climbing-image NEB with a rotating-orbit refinement for
  locating complex transition states.
* **Machine-Learned Potentials with `metatomic`**: Full support for `metatomic`
  potentials via the [metatensor
    ecosystem](https://docs.metatensor.org/latest/index.html).
* **Expanded Potential Interfaces**: Drastically increased interoperability with
  direct interfaces for **NWChem**, **ORCA**, **AMS**, **XTB**, and **ASE**.
* **Input Validation and Schema**: A Pydantic schema for all
  configuration files provides automatic input validation and clearer error
  messages.
* **Meson build**: one build on Linux, macOS, and Windows.

Papers that use these methods are listed with the rest of the changes:

```{toctree}
:maxdepth: 2
:caption: Contents

release-notes
publications
```
