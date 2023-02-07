---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.4
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

```{code-cell} ipython3
:tags: []

import numpy as np
import numpy.typing as npt

import pyeonclient as ec
```

```{code-cell} ipython3
:tags: []

params = ec.Parameters()
params.load("../../examples/neb-al/config.ini") # NEB-Al
ec.log_init(params, "blah.log") # Just in case
```

```{code-cell} ipython3
:tags: []

product = ec.Matter(params)
product.con2matter("../../examples/neb-al/product.con")
reactant = ec.Matter(params)
reactant.con2matter("../../examples/neb-al/reactant.con")
```

```{code-cell} ipython3
:tags: []

product.forces
```

```{code-cell} ipython3
:tags: []

product.pot_energy
```

```{code-cell} ipython3
:tags: []

pot = ec.Morse()
```

```{code-cell} ipython3
:tags: []

type(pot.getPotential(params))
```

```{code-cell} ipython3
## Broken sometimes
pot.force(product)
```

```{code-cell} ipython3
:tags: []

pot.force(product.numberOfAtoms(), product.positions, np.ones(product.numberOfAtoms()), product.getCell())
```

```{code-cell} ipython3
:tags: []

## Works in the terminal, and ipython, but crashes here..
#lldb python -- -c 'import pyeonclient as ec; params = ec.Parameters(); params.load("../examples/neb-al/config.ini"); ec.log_init(params, "blah.log"); m1 = ec.Matter(params); m1.con2matter("../examples/neb-al/product.con"); pot = ec.Morse(); print(m1.pot_energy); m1.setPotential(pot); print(m1.pot_energy)'
## Kernel restart here
#product.setPotential(pot)
#product.pot_energy
```

```{code-cell} ipython3
:tags: []

m1 = ec.Matter(params); m1.con2matter("../../examples/neb-al/product.con"); pot.force(m1);
```

```{code-cell} ipython3
:tags: []

type(pot.getPotential(params))
```

```{code-cell} ipython3
pot.force(m1)
```

```{code-cell} ipython3
:tags: []

## Broken
m1.setPotential(pot.getPotential(params))
m1.pot_energy
```

```{code-cell} ipython3
:tags: []

class FakePot(ec.Potential):
    def force(self, a, b, c, d, e, f, g):
        return 32
    def initialize(self):
        pass
```

```{code-cell} ipython3
:tags: []

fp = FakePot()
```

```{code-cell} ipython3
:tags: []

m1.setPotential(fp)
```

```{code-cell} ipython3
:tags: []

m1.pot_energy
```

```{code-cell} ipython3
:tags: []

m1.getPotential()
```

```{code-cell} ipython3

```
