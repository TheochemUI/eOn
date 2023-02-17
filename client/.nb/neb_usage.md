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

pot = ec.makePotential(params)
```

```{code-cell} ipython3
:tags: []

product.pot_energy
```

```{code-cell} ipython3
:tags: []

neb = ec.NudgedElasticBand(reactant, product, params)
neb.printImageData(True)
```

```{code-cell} ipython3
:tags: []

neb.compute()
```

```{code-cell} ipython3
:tags: []

nmatter = neb.neb_images
```

```{code-cell} ipython3
:tags: []

r2 = ec.Matter(params)
r2.con2matter("../../examples/neb-al/reactant.con")
r2.positions += 0.005
```

```{code-cell} ipython3
:tags: []

nmatter[3] = r2
```

```{code-cell} ipython3
:tags: []

neb.neb_images = nmatter
```

```{code-cell} ipython3
:tags: []

neb.compute()
```

```{code-cell} ipython3
:tags: []

neb.maxEnergyImage
```

```{code-cell} ipython3
:tags: []

neb.extremumCurvature
```

```{code-cell} ipython3
dir(neb)
```
