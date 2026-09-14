---
myst:
  html_meta:
    "description": "Public C++ library headers for embedding eOn NEB, IDPP, and dimer."
    "keywords": "eOn library API, gpr_optim, NEB, IDPP, dimer"
---

# Library API

External C++ (gpr_optim, a custom driver) should include `eon/api.h`.
That header pulls Matter, Parameters, Potential, NEB, IDPP/SIDPP path
builders, ImprovedDimer, Lanczos, and JobResult. Job implementations
(`*Job.h`) stay internal.

```cpp
#include "eon/api.h"

eonc::Parameters params;
params.potential_options.potential = eonc::PotType::LJ;
auto pot = eonc::helpers::makePotential(params);
auto path = eonc::helpers::neb_paths::sidppPath(reactant, product, 5, params);
eonc::NudgedElasticBand neb(path, params, pot);
neb.compute();
```
