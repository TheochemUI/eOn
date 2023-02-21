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
```

```{code-cell} ipython3
:tags: []

m1 = ec.Matter(params)
```

```{code-cell} ipython3
:tags: []

m1.con2matter("../gtests/data/systems/one_Pt_on_frozenSurface/pos.con")
```

```{code-cell} ipython3
:tags: []

m1.pot_energy
```

```{code-cell} ipython3
m1.getFreeV()
```

```{code-cell} ipython3
:tags: []

m1.getFree()
```

```{code-cell} ipython3
:tags: []

m1.positions[0]
```

```{code-cell} ipython3
:tags: []

mpot = ec.Morse(params)
```

```{code-cell} ipython3
:tags: []

mpot.get_ef(m1.positions, np.ones(m1.numberOfAtoms()), m1.getCell())
```

For the Pt system with one moving atom, we will plot the section bound by $X \in [7.6273, 11.66]$ and $Y \in [7.63402, 13.074]$ and a suitably constrained Z value, $14.584501$ in this case.

+++ {"tags": []}

We will reuse plotting components from `anneal`.

```{code-cell} ipython3
:tags: []

# TODO
class CorrectedForces:
    """
    This class is meant to take a matter object as its constructor and correct for masks (moving, frozen).
    It is the callers responsibility to ensure that the function is
    called on either a matter object with the same number of atoms or
    that the positions are for free atoms only.
    """
    def __init__(self, ec_pot: ec.Potential, ec_matter: ec.Matter):
        self.pot = ec_pot
        self.mask = pos[np.mean(self.freeIdx, axis=1, dtype=bool)] 
    def operator():
        pass
```

```{code-cell} ipython3
:tags: []

from anneal.core.components import ObjectiveFunction, NumLimit, FPair
from anneal.viz.viz2d import Plot2dObj
import numpy as np
```

```{code-cell} ipython3
:tags: []

class EONPot1Pt(ObjectiveFunction):
    """Just constrained to dealing with 1 moving atom constrained in 2D
    It is assumed here, that the positions are only free positions"""
    def __init__(
        self, ec_pot, ec_mat, limits=NumLimit(dims=2,
                              low=np.array([7.6273, 7.63402, 14.584501]),
                              high=np.array([11.66, 13.074, 14.584501]))
    ):
        self.pot = ec_pot
        self.mat = ec_mat
        self.freeIdx = ec_mat.getFree()
        self.freeMask = np.mean(self.freeIdx, axis=1, dtype=bool)
        super().__init__(
            limits,
            None # Don't know the minimum
        )

    def singlepoint(self, pos):
        """Here this is a 'matter style' single point"""
        cpos = self.mat.positions
        effective_pos= np.append(pos.ravel(), 14.584501)
        cpos[self.freeMask] = effective_pos
        self.mat.positions = cpos
        return self.mat.pot_energy

    def multipoint(self, pos: list):
        """Returns a corrected matter approximation, not suitable for multiple matter objects"""
        elist = []
        self.mat.setPotential(self.pot)
        elist = np.apply_along_axis(self.singlepoint, axis=1, arr=pos.reshape(-1, 2))
        return elist

    def __repr__(self):
        return "EON 1 free atom Holder"
```

```{code-cell} ipython3
:tags: []

epp = EONPot1Pt(mpot, m1)
```

```{code-cell} ipython3
:tags: []

epp(np.array([1,2])); print(epp(np.arange(90).reshape(-1, 2)));
```

```{code-cell} ipython3
:tags: []

eplot = Plot2dObj(epp, 30)
```

```{code-cell} ipython3
:tags: []

eplot.createContour(showGlob=False)
```

```{code-cell} ipython3
:tags: []

eplot.create3d(showGlob=False)
```

```{code-cell} ipython3
:tags: []

from anneal.quenchers.boltzmann import BoltzmannQuencher
from anneal.core.components import AcceptStates, Quencher
```

```{code-cell} ipython3
:tags: []

##bq = BoltzmannQuencher(epp, T_init=5)
```

```{code-cell} ipython3

```
