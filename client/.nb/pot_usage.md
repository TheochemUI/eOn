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

mpot = ec.Morse(params)
```

```{code-cell} ipython3
:tags: []

product.setPotential(mpot)
```

```{code-cell} ipython3
:tags: []

product.pot_energy
```

```{code-cell} ipython3
#mpot.get_ef(product)
```

```{code-cell} ipython3
:tags: []

mpot.get_ef(product.positions, np.ones(product.numberOfAtoms()), product.getCell())
```

```{code-cell} ipython3
:tags: []

m1 = ec.Matter(params); m1.con2matter("../../examples/neb-al/product.con"); #mpot.force(m1);
```

```{code-cell} ipython3
#mpot.force(m1)[1].shape
```

```{code-cell} ipython3
:tags: []

class FakePot(ec.Potential):
    def __init__(self, params):
        #params.potential = ec.PotType.PYTHON
        ec.Potential.__init__(self, params)
    def get_ef(self, a, b, c):
        return (33, np.ones(601*3).reshape(601, 3))
```

```{code-cell} ipython3
:tags: []

fp = FakePot(params)
```

```{code-cell} ipython3
:tags: []

m1.setPotential(fp)
```

```{code-cell} ipython3
:tags: []

fp.force(m1)
```

```{code-cell} ipython3
:tags: []

m1.pot_energy
```

```{code-cell} ipython3
:tags: []

m1.forces
```

```{code-cell} ipython3
:tags: []

##ec.callPotential(fp, product.positions, np.ones(product.numberOfAtoms()), product.getCell())
```

```{code-cell} ipython3
:tags: []

params.potential = ec.PotType.LJ
ljpot = ec.makePotential(params)
ljpot.get_ef(product.positions, np.ones(product.numberOfAtoms()), product.getCell())
```
