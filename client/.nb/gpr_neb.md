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

import torch
import gpytorch
import math
import copy
from matplotlib import cm
from matplotlib import pyplot as plt
import numpy as np

from scipy import spatial

%matplotlib inline
#%load_ext autoreload
#%autoreload 2
```

```{code-cell} ipython3
:tags: []

import pyeonclient as ec
import eoncpy as ecp
```

```{code-cell} ipython3
:tags: []

from eindir.core.components import *
from eindir.viz.viz2d import *
```

# GPR Implmentations for EON
This is the final class based structure for re-training and using GPRs.

```{code-cell} ipython3
:tags: []

def make_masks(ecmatter, radius=3.5):
    freeIdx = ecmatter.getFree()
    freeMask = np.mean(freeIdx, axis=1, dtype=bool)
    tree1 = spatial.KDTree(ecmatter.positions)
    tree2 = spatial.KDTree(ecmatter.free_positions)
    freeMoveMask = copy.deepcopy(freeMask)
    for neigh_idx in tree2.query_ball_tree(tree1, r=radius):
        freeMoveMask[neigh_idx] = True
    return (freeMask, freeMoveMask)
```

```{code-cell} ipython3
:tags: []

def make_data_free_move(listMatter, radius = 3.5):
    freeMask, freeMoveMask = make_masks(listMatter[0], radius)
    listPos = [x.positions[freeMoveMask] for x in listMatter]
    listForce = [x.forces[freeMoveMask] for x in listMatter]
    listEnergy = np.array([x.pot_energy for x in listMatter])
    x_train = torch.from_numpy(
        np.row_stack([x.ravel() for x in listPos])
    )
    y_train = torch.from_numpy(
        np.row_stack([np.append(x, y) for x,y in zip(listEnergy, listForce)])
    )
    return (x_train, y_train)
```

## Usage
In practice, this should be straightforward.

```{code-cell} ipython3
:tags: []

params = ec.Parameters()
p2 = ec.Parameters()
params.load("../../client/tests/neb/config.ini") # Morse Pt
p2.load("../../client/tests/neb/config.ini")
ec.log_init(params, "blah.log") # Just in case
product = ec.Matter(params)
product.con2matter("../../client/tests/neb/product.con")
reactant = ec.Matter(params)
reactant.con2matter("../../client/tests/neb/reactant.con")
m1 = ec.Matter(params)
m1.con2matter("../../client/tests/neb/reactant.con")
```

```{code-cell} ipython3
:tags: []

neb = ec.NudgedElasticBand(reactant, product, params)
```

```{code-cell} ipython3
:tags: []

initData = neb.neb_images
xtrain, ytrain = make_data_free_move(initData)
```

```{code-cell} ipython3
:tags: []

mask = ecp.structures.MatterMask(reactant)
mgpr = ecp.GPR.MakeGPR(params, initData, mask, 0)
```

```{code-cell} ipython3
:tags: []

mgpr.train(70, log=True)
```

```{code-cell} ipython3
:tags: []

epot = mgpr.get_ecpot(p2)
#ecp.potentials.test_force_convergence(params, mgpr.masks, neb.neb_images, gneb.mgpr.get_ecpot(p2))[1:-1]
ecp.potentials.test_force_convergence(params, mgpr.masks, neb.neb_images, epot)[1:-1]
```

```{code-cell} ipython3
:tags: []

#neb.compute()
```

```{code-cell} ipython3
:tags: []

ec.makePotential(params)
```

```{code-cell} ipython3
:tags: []

mpot = ec.Morse(params)
```

```{code-cell} ipython3
:tags: []

epot.get_ef(product.positions, m1.numberOfAtoms, m1.getCell())
```

```{code-cell} ipython3
:tags: []

[x.pot_energy for x in initData]
```

```{code-cell} ipython3
:tags: []

class GPRNeb:
    def __init__(
        self,
        gprmk: ecp.GPR.MakeGPR,
        reactant: ec.Matter,
        product: ec.Matter,
        params: ec.Parameters,
    ):
        self.mgpr = gprmk
        self.reactant = reactant
        self.product = product
        self.mat_params = params
        self.nrounds = 0
        self.conv_rounds = []
        self.gprpot = gprmk.get_ecpot(ec.Parameters())

    def run_neb_round(self):
        neb = ec.NudgedElasticBand(self.reactant, self.product, self.mat_params)
        [x.setPotential(self.gprpot) for x in neb.neb_images]
        neb.compute()
        conv_data = ecp.potentials.test_force_convergence(self.mat_params,
                                                     self.mgpr.masks,
                                                     neb.neb_images, self.gprpot)[1:-1]
        if max(conv_data) <= 0.5:
            return
        else:
            self.mgpr.add_observations(neb.neb_images)
            self.gprpot = self.mgpr.get_ecpot(ec.Parameters())
            self.nrounds += 1
            self.conv_rounds.append(conv_data)
        return
```

```{code-cell} ipython3
:tags: []

gneb = GPRNeb(mgpr, reactant, product, params)
```

```{code-cell} ipython3
:tags: []

%timeit
gneb.run_neb_round()
```

```{code-cell} ipython3
:tags: []

gneb.conv_rounds
```

```{code-cell} ipython3
:tags: []

class EONPot1Pt(ObjectiveFunction):
    """Just constrained to dealing with 1 moving atom constrained in 2D
    It is assumed here, that the positions are only free positions.
    Additionally, the Z axis is constrained, so the dimensions are correct for SA"""
    def __init__(
        self, ec_pot, ec_mat, limits=NumLimit(dims=2,
                              low=np.array([5.6273, 5.63402]),#, 14.584501]),
                              high=np.array([11.66, 13.074]))#, 14.584501]))
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
        effective_pos= np.append(pos.ravel(), 14.573576)
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

gneb.nrounds
```

```{code-cell} ipython3
:tags: []

epp = EONPot1Pt(gneb.mgpr.get_ecpot(p2), m1)
eplt = Plot2dObj(epp, 40)
eplt.createContour()
```

```{code-cell} ipython3
:tags: []

epp = EONPot1Pt(mpot, m1)
eplt = Plot2dObj(epp, 40)
eplt.createContour()
```

```{code-cell} ipython3
:tags: []

neb.extremumCurvature
```

```{code-cell} ipython3
:tags: []

neb.compute()
```

```{code-cell} ipython3
:tags: []

x_train, y_train = make_data_free_move(neb.neb_images)
```

```{code-cell} ipython3
:tags: []

mgpr.add_observations(neb.neb_images)
```

```{code-cell} ipython3
:tags: []

epot2 = mgpr.get_ecpot(p2)
```

```{code-cell} ipython3
:tags: []

neb2 = ec.NudgedElasticBand(reactant, product, params)
```

```{code-cell} ipython3
:tags: []

[x.setPotential(epot2) for x in neb2.neb_images]
```

```{code-cell} ipython3
:tags: []

neb2.compute()
```

```{code-cell} ipython3
:tags: []

ecp.potentials.test_force_convergence(params, mgpr.masks, neb.neb_images, epot)
```

```{code-cell} ipython3
:tags: []

ecp.potentials.test_force_convergence(params, mgpr.masks, neb2.neb_images, epot2)
```

```{code-cell} ipython3
:tags: []

[x.pot_energy for x in neb2.neb_images]
```

```{code-cell} ipython3
:tags: []

mgpr.add_observations(neb2.neb_images)
```

```{code-cell} ipython3
:tags: []

epot3 = mgpr.get_ecpot(p2)
```

```{code-cell} ipython3
:tags: []

ecp.potentials.test_force_convergence(params, mgpr.masks, neb2.neb_images, epot3)
```

```{code-cell} ipython3
:tags: []

neb3 = ec.NudgedElasticBand(reactant, product, params)
[x.setPotential(epot3) for x in neb3.neb_images]
```

```{code-cell} ipython3
:tags: []

neb3.compute()
```

```{code-cell} ipython3
:tags: []

ecp.potentials.test_force_convergence(params, mgpr.masks, neb3.neb_images, epot3)
```

```{code-cell} ipython3
neb3.extremumCurvature
```

```{code-cell} ipython3
:tags: []

neb_true = ec.NudgedElasticBand(reactant, product, params)
neb_true.compute()
```

```{code-cell} ipython3
:tags: []

ecp.potentials.test_force_convergence(params, mgpr.masks, neb_true.neb_images, epot3)
```

```{code-cell} ipython3
:tags: []

ecp.potentials.test_force_convergence(params, mgpr.masks, neb3.neb_images, epot3)
```

```{code-cell} ipython3
:tags: []

[x for x in mgpr.model.train_targets]
```

```{code-cell} ipython3
:tags: []

freeMask, freeMoveMask = make_masks(reactant, radius=3.5)
```

```{code-cell} ipython3
:tags: []

from eoncpy import structures
```

```{code-cell} ipython3
:tags: []

mm = structures.Masks(reactant.positions, np.mean(reactant.getFree(), axis=1).astype(bool))
```

```{code-cell} ipython3
:tags: []

mm.mark_neighbors_of_free(3.5)
```

```{code-cell} ipython3
:tags: []

np.all(mm.is_active == freeMoveMask)
```

```{code-cell} ipython3
:tags: []

pot = GPRPotFreeMove(params, gp1.model, gp1.likelihood,
                     copy.deepcopy(freeMoveMask),
                     copy.deepcopy(freeMask))
```

```{code-cell} ipython3
:tags: []

for img in nebgp.neb_images:
    img.setPotential(pot)
    print(img.pot_energy)
```

```{code-cell} ipython3
:tags: []

reactant.setPotential(pot)
product.setPotential(pot)
nebgp = ec.NudgedElasticBand(reactant, product, params)
```

```{code-cell} ipython3
:tags: []

[x.free_positions for x in nebgp.neb_images]
```

```{code-cell} ipython3
:tags: []

xtrain, ytrain = get_more_data(xtrain, ytrain, nextData)
```

```{code-cell} ipython3
:tags: []

gp2 = makeGPR(xtrain, ytrain)
```

```{code-cell} ipython3
:tags: []

gp2.train(log=False)
```

```{code-cell} ipython3
:tags: []

pot2 = GPRPotFreeMove(params, gp2.model, gp2.likelihood,
                     copy.deepcopy(freeMoveMask),
                     copy.deepcopy(freeMask))
```

```{code-cell} ipython3
:tags: []

epp = EONPot1Pt(epot, m1)
eplot = Plot2dObj(epp, 30)
eplot.createContour(showGlob=True)
```

```{code-cell} ipython3
:tags: []

m1.setPotential(pot)
m1.pot_energy
```

```{code-cell} ipython3
:tags: []

m1.setPotential(pot2)
m1.pot_energy
```

```{code-cell} ipython3
:tags: []

 get_more_data(xtrain, ytrain, nextData)[0].shape
```

```{code-cell} ipython3
:tags: []

[x.pot_energy for x in nextData]
```

```{code-cell} ipython3

```
