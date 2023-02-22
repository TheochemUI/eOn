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
from matplotlib import cm
from matplotlib import pyplot as plt
import numpy as np

from scipy.spatial import KDTree

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


We will first setup the data for this. As a first approximation, we will train a new GPR at each stage, later we can iterate on this design. Taking all the atoms into consideration is too expensive, so only free atoms are considered.

```{code-cell} ipython3
:tags: []

initData = neb.neb_images
```

```{code-cell} ipython3
:tags: []

def make_data_free(listMatter):
    listPos = [x.free_positions for x in neb.neb_images]
    listForce = [x.free_forces for x in neb.neb_images]
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

# x_train = torch.from_numpy(
#     np.row_stack([
#         Xpos_init[0].ravel(),
#         Xpos_init[1].ravel(),]))   
# y_train = torch.from_numpy(
#     np.row_stack([
#     np.append(Ye_init_np[0], Xfrcs_init[0].ravel()),
#     np.append(Ye_init_np[1], Xfrcs_init[1].ravel()),
# ]))## should be 1, ndim+1
x_train, y_train = make_data_free(initData)
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

training_iter = 200
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

```{code-cell} ipython3
:tags: []

# # Set into eval mode
# model.eval()
# likelihood.eval()
# # Predict on a new point
# with torch.no_grad(), gpytorch.settings.fast_computations(log_prob=False, covar_root_decomposition=False):
#     test_x = torch.from_numpy(Xpos_init[2].ravel()).unsqueeze(0)
#     predictions = likelihood(model(test_x))
#     mean = predictions.mean
```

```{code-cell} ipython3
:tags: []

# frcs = mean.squeeze(0)[1::1].reshape(-1, 3)
# energy = mean.squeeze(0)[0]
```

## Setup Single Matter <-> GPR Model
The rational behind this is in `gpr_eon.md`.

```{code-cell} ipython3
:tags: []

class GPRPot(ec.Potential):
    def __init__(self, params, model, likelihood):
        self.model = model
        self.likelihood = likelihood
        ec.Potential.__init__(self, params)
    def get_ef(self, pos, atmnrs, box):
        self.model.eval()
        self.likelihood.eval()
        with torch.no_grad(), gpytorch.settings.fast_computations(log_prob=False,
                                                                 covar_root_decomposition=False):
            test_x = torch.from_numpy(pos.ravel()).unsqueeze(0)
            predictions = self.likelihood(self.model(test_x))
            mean = predictions.mean
        frcs = mean.squeeze(0)[1::1].reshape(-1, 3)
        energy = mean.squeeze(0)[0]
        return (energy, frcs)
```

```{code-cell} ipython3
:tags: []

class GPRPotFree(ec.Potential):
    def __init__(self, params, model, likelihood, freeIdx):
        self.model = model
        self.freeIdx = freeIdx
        self.likelihood = likelihood
        ec.Potential.__init__(self, params)
    def get_ef(self, pos, atmnrs, box):
        self.model.eval()
        self.likelihood.eval()
        free_pos = pos[np.mean(self.freeIdx, axis=1, dtype=bool)]
        with torch.no_grad(), gpytorch.settings.fast_computations(log_prob=False,
                                                                 covar_root_decomposition=False):
            test_x = torch.from_numpy(free_pos.ravel()).unsqueeze(0)
            predictions = self.likelihood(self.model(test_x))
            mean = predictions.mean
        frcs = mean.squeeze(0)[1::1].reshape(-1, 3)
        ## Need to set the fixed forces to zero
        forces = np.zeros_like(pos)
        forces[np.mean(self.freeIdx, axis=1, dtype=bool)] = frcs
        ## free_forces = frcs[np.mean(self.freeIdx, axis=1, dtype=bool)]
        energy = mean.squeeze(0)[0]
        return (energy, forces)
```

```{code-cell} ipython3
:tags: []

pot = GPRPotFree(params, model, likelihood, reactant.getFree())
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

datTwo = initData + neb_gpr.neb_images
```

```{code-cell} ipython3
:tags: []

x_train, y_train = make_data_free(datTwo)
```

```{code-cell} ipython3
:tags: []

neb.compute()
```

```{code-cell} ipython3
:tags: []

neb.compute()
```

```{code-cell} ipython3
:tags: []

neb_gpr.extremumPosition
```

```{code-cell} ipython3
:tags: []

neb.extremumPosition
```

```{code-cell} ipython3
:tags: []

neb.extremumCurvature
```

```{code-cell} ipython3
:tags: []

neb_gpr.extremumCurvature
```

```{code-cell} ipython3
:tags: []

reactant.setPotential(pot)
```

```{code-cell} ipython3
:tags: []

reactant.pot_energy
```

```{code-cell} ipython3
:tags: []

reactant.forces
```

```{code-cell} ipython3
:tags: []

[x.pot_energy for x in neb.neb_images]
```

```{code-cell} ipython3

```
