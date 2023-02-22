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
%load_ext autoreload
%autoreload 2
```

```{code-cell} ipython3
:tags: []

import pyeonclient as ec
```

```{code-cell} ipython3
:tags: []

params = ec.Parameters()
params.load("../../client/tests/neb/config.ini") # Morse Pt
ec.log_init(params, "blah.log") # Just in case
```

```{code-cell} ipython3
:tags: []

product = ec.Matter(params)
product.con2matter("../../client/tests/neb/product.con")
reactant = ec.Matter(params)
reactant.con2matter("../../client/tests/neb/reactant.con")
```

```{code-cell} ipython3
:tags: []

neb = ec.NudgedElasticBand(reactant, product, params)
```

# GPR Potential

This is a first pass at making the standard GPR model with a squared exponential.

+++

## Setup Data
We will first setup the data for this. As a first approximation, we will train a new GPR at each stage, later we can iterate on this design. Taking all the atoms into consideration is too expensive, so only free atoms were considered originally. However, the resulting potential is poor. This means we must re-create the moving / frozen aspects of the original formulation.

We will use the `kdtree` from `scipy`.

```{code-cell} ipython3
:tags: []

initData = neb.neb_images
```

```{code-cell} ipython3
:tags: []

m1 = reactant
```

```{code-cell} ipython3
:tags: []

freeIdx = m1.getFree()
freeMask = np.mean(freeIdx, axis=1, dtype=bool)
```

```{code-cell} ipython3
:tags: []

m1.numberOfFreeAtoms == np.apply_along_axis(np.mean, axis=1, arr=m1.getFree().reshape(-1, 3)).sum()
```

```{code-cell} ipython3
:tags: []

tree1 = spatial.KDTree(m1.positions)
```

```{code-cell} ipython3
:tags: []

tree2 = spatial.KDTree(m1.free_positions)
```

```{code-cell} ipython3
:tags: []

tree2.query_ball_tree(tree1, r=3.5)
```

So now we need to set a mask on the positions which includes the free atoms, along with these, the neighbors within 3.5 units of a moving atom.

```{code-cell} ipython3
:tags: []

freeMoveMask = copy.deepcopy(freeMask)
for neigh_idx in tree2.query_ball_tree(tree1, r=3.5):
    freeMoveMask[neigh_idx] = True
```

```{code-cell} ipython3
:tags: []

m1.positions[freeMoveMask]
```

```{code-cell} ipython3
:tags: []

def make_data_free_move(listMatter):
    listPos = [x.positions[freeMoveMask] for x in neb.neb_images]
    listForce = [x.forces[freeMoveMask] for x in neb.neb_images]
    listEnergy = np.array([x.pot_energy for x in neb.neb_images])
    x_train = torch.from_numpy(
        np.row_stack([x.ravel() for x in listPos])
    )
    y_train = torch.from_numpy(
        np.row_stack([np.append(x, y) for x,y in zip(listEnergy, listForce)])
    )
    return (x_train, y_train)
```

```{code-cell} ipython3
:tags: []

x_train, y_train = make_data_free_move(initData)
```

```{code-cell} ipython3
:tags: []

y_train.shape, x_train.shape
```

```{code-cell} ipython3
:tags: []

class GPModelWithDerivativesTwo(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood): 
        super(GPModelWithDerivativesTwo, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMeanGrad()
        self.base_kernel = gpytorch.kernels.RBFKernelGrad(ard_num_dims=train_x.shape[1]) ## Important, additional dimensions
        self.covar_module = gpytorch.kernels.ScaleKernel(self.base_kernel)

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultitaskMultivariateNormal(mean_x, covar_x)
```

```{code-cell} ipython3
:tags: []

likelihood = gpytorch.likelihoods.MultitaskGaussianLikelihood(num_tasks=y_train.shape[1]) # Value + x-derivative + y-derivative + z-derivative
model = GPModelWithDerivativesTwo(x_train,
                               y_train,
                                 likelihood)
```

```{code-cell} ipython3
:tags: []

training_iter = 400
# Find optimal model hyperparameters
model.train()
likelihood.train()

# Use the adam optimizer
optimizer = torch.optim.Adam(model.parameters(), lr=0.05)  # Includes GaussianLikelihood parameters

# "Loss" for GPs - the marginal log likelihood
mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)

for i in range(training_iter):
    optimizer.zero_grad()
    output = model(x_train)
    loss = -mll(output, y_train)
    loss.backward()
    print("Iter %d/%d - Loss: %.3f   lengthscales: %.3f, %.3f   noise: %.3f" % (
        i + 1, training_iter, loss.item(),
        model.covar_module.base_kernel.lengthscale.squeeze()[0],
        model.covar_module.base_kernel.lengthscale.squeeze()[1],
        model.likelihood.noise.item()
    ))
    optimizer.step()
```

## Setup Single Matter <-> GPR Model
The rational behind this is in `gpr_eon.md`.

```{code-cell} ipython3
:tags: []

class GPRPotFreeMove(ec.Potential):
    def __init__(self, params, model, likelihood, move_mask, free_mask):
        self.model = model
        self.move_mask = move_mask
        self.free_mask = free_mask
        self.likelihood = likelihood
        ec.Potential.__init__(self, params)
    def get_ef(self, pos, atmnrs, box):
        self.model.eval()
        self.likelihood.eval()
        free_move_pos = pos[self.move_mask]
        with torch.no_grad(), gpytorch.settings.fast_computations(log_prob=False,
                                                                 covar_root_decomposition=False):
            test_x = torch.from_numpy(free_move_pos.ravel()).unsqueeze(0)
            predictions = self.likelihood(self.model(test_x))
            mean = predictions.mean
        frcs = mean.squeeze(0)[1::1].reshape(-1, 3)
        ## Need to set the fixed forces to zero
        forces = np.zeros_like(pos)
        forces[self.move_mask] = frcs
        forces[self.free_mask] = 0
        energy = mean.squeeze(0)[0]
        return (energy, forces)
```

```{code-cell} ipython3
:tags: []

pot = GPRPotFreeMove(params, model, likelihood, copy.deepcopy(freeMoveMask), copy.deepcopy(freeMask))
```

```{code-cell} ipython3
:tags: []

neb_gpr = ec.NudgedElasticBand(reactant, product, params, pot)
neb_gpr.printImageData(True)
```

```{code-cell} ipython3
:tags: []

neb_gpr.compute()
```

```{code-cell} ipython3
:tags: []

from anneal.core.components import ObjectiveFunction, NumLimit, FPair, MAX_LIMITS
from anneal.viz.viz2d import Plot2dObj

class EONPot1Pt(ObjectiveFunction):
    """Just constrained to dealing with 1 moving atom constrained in 2D
    It is assumed here, that the positions are only free positions.
    Additionally, the Z axis is constrained, so the dimensions are correct for SA"""
    def __init__(
        self, ec_pot, ec_mat, limits=NumLimit(dims=2,
                              low=np.array([-20, -20]),#, 14.584501]),
                              high=np.array([20, 20]))#, 14.584501]))
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
        effective_pos= np.append(pos.ravel(), 12.296689)
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

epp = EONPot1Pt(pot, m1)
```

```{code-cell} ipython3
:tags: []

eplot = Plot2dObj(epp, 30)
```

```{code-cell} ipython3
:tags: []

eplot.createContour(showGlob=True)
```

```{code-cell} ipython3
:tags: []

pot.get_ef(product.positions, m1.numberOfAtoms, m1.getCell())
```

```{code-cell} ipython3
:tags: []

product.pot_energy
```

It seems that the only difference is in the forces, but the difference is staggering. Much better to have the moving atoms mask.

```{code-cell} ipython3

```
