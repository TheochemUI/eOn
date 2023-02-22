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

neb = ec.NudgedElasticBand(reactant, product, params)
```

```{code-cell} ipython3
:tags: []

Xpos_init = [x.positions for x in neb.neb_images]
Xfrcs_init = [x.forces for x in neb.neb_images]
Yenerg_init = [x.pot_energy for x in neb.neb_images]
```

```{code-cell} ipython3
:tags: []

neb.compute()
```

```{code-cell} ipython3
:tags: []

Xpos_fin = [x.positions for x in neb.neb_images]
Xfrcs_fin = [x.positions for x in neb.neb_images]
Yenerg_fin = [x.pot_energy for x in neb.neb_images]
```

```{code-cell} ipython3
:tags: []

## 1D variants
Xpos_rinit = [x.ravel() for x in Xpos_init]
Xfrcs_rinit = [x.ravel() for x in Xfrcs_init]
Xpos_rfin = [x.ravel() for x in Xpos_fin]
Xfrcs_rfin = [x.ravel() for x in Xfrcs_fin]
```

```{code-cell} ipython3
:tags: []

## NumPy variants
Xp_init_np = np.concatenate([Xpos_rinit], axis=1)
Xf_init_np = np.concatenate([Xfrcs_rinit], axis=1)
Xp_fin_np = np.concatenate([Xpos_rfin], axis=1)
Xf_fin_np = np.concatenate([Xfrcs_rfin], axis=1)
Ye_init_np = np.array(Yenerg_init)
Ye_fin_np = np.array(Yenerg_fin)
```

```{code-cell} ipython3
:tags: []

## PyTorch variants
Xpi = torch.from_numpy(Xp_init_np)
Xfi = torch.from_numpy(Xf_init_np)
Xpf = torch.from_numpy(Xp_fin_np)
Xff = torch.from_numpy(Xf_fin_np)
Yei = torch.from_numpy(Ye_init_np)
Yef = torch.from_numpy(Ye_fin_np)
```

```{code-cell} ipython3
:tags: []

to_torchvec = lambda x: torch.from_numpy(np.concatenate([x], axis=1))
```

```{code-cell} ipython3
:tags: []

np.concatenate([Xpos_init], axis=1).shape
```

## Single Matter to multiple GPR data
Essentially we treat the system as a collection of atomic data records, that is one matter image provides us with N observations, each of the kind (Energy, Dx, Dy).

Importantly, this means that the GPR only outputs 4 terms, the 1D tuple (E, dx, dy, dz), and that the input dimensionality of the GPR is essentially 3, with the likelihood having 4 terms.

```{code-cell} ipython3
:tags: []

train_x = torch.tensor(Xpos_init)
train_xf = torch.tensor(Xfrcs_init)
```

```{code-cell} ipython3
:tags: []

x_train = train_x[0]
xf_train = train_xf[0]
train_y = torch.ones(601) * Yei[0]
```

```{code-cell} ipython3
:tags: []

train_y.shape, x_train[:,1].shape
```

```{code-cell} ipython3
:tags: []

y_train = torch.stack([train_y, x_train[:, 0], x_train[:, 1], x_train[:, 2]], -1).squeeze(1)
```

```{code-cell} ipython3
:tags: []

x_train.shape, y_train.shape
```

```{code-cell} ipython3
:tags: []

ytt = torch.stack([torch.ones(601)*Yei[1], train_x[1][:, 0], train_x[1][:, 1], train_x[2][:, 2]], -1).squeeze(1)
```

```{code-cell} ipython3
:tags: []

ytrain=torch.vstack((y_train, ytt))
```

```{code-cell} ipython3
:tags: []

xtrain = torch.vstack((train_x[0], train_x[1]))
```

```{code-cell} ipython3
:tags: []

class GPModelWithDerivatives(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood):
        super(GPModelWithDerivatives, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMeanGrad()
        self.base_kernel = gpytorch.kernels.RBFKernelGrad(ard_num_dims=3) ## Important, additional dimensions
        self.covar_module = gpytorch.kernels.ScaleKernel(self.base_kernel)

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultitaskMultivariateNormal(mean_x, covar_x)

likelihood = gpytorch.likelihoods.MultitaskGaussianLikelihood(num_tasks=4)  # Value + x-derivative + y-derivative + z-derivative
```

```{code-cell} ipython3
:tags: []

model = GPModelWithDerivatives(xtrain,
                               ytrain,
                               likelihood)
```

```{code-cell} ipython3
:tags: []

training_iter = 50
# Find optimal model hyperparameters
model.train()
likelihood.train()

# Use the adam optimizer
optimizer = torch.optim.Adam(model.parameters(), lr=0.05)  # Includes GaussianLikelihood parameters

# "Loss" for GPs - the marginal log likelihood
mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)

for i in range(training_iter):
    optimizer.zero_grad()
    output = model(xtrain)
    loss = -mll(output, ytrain)
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

# Set into eval mode
model.eval()
likelihood.eval()
```

```{code-cell} ipython3
:tags: []

train_x[0]
```

```{code-cell} ipython3
:tags: []

with torch.no_grad(), gpytorch.settings.fast_computations(log_prob=False, covar_root_decomposition=False):
    test_x = train_x[2]
    predictions = likelihood(model(test_x))
    mean = predictions.mean
```

```{code-cell} ipython3
:tags: []

predictions
```

```{code-cell} ipython3
:tags: []

mean[:, 0].mean()
```

```{code-cell} ipython3
:tags: []

Xfrcs_init[1]
```

```{code-cell} ipython3
:tags: []

def abs_diff_terms(exp, true_frcs, true_e):
    return (
        np.abs(exp[:, 0]) - np.abs(true_e),
        np.abs(exp[:, 1]) - np.abs(true_frcs[:, 0]),
        np.abs(exp[:, 2]) - np.abs(true_frcs[:, 1]),
        np.abs(exp[:, 3]) - np.abs(true_frcs[:, 2]),
    )
```

```{code-cell} ipython3
:tags: []

fe, fdx, fdy, fdz = abs_diff_terms(mean, Xfrcs_init[2], Yenerg_init[2])
print(f"Absolute difference in forces\n: x = {fdx.mean()}\n y = {fdy.mean()}\n z = {fdz.mean()}\n the energy difference is {fe.mean()}")
```

## Single Matter to single GPR data
Essentially we treat the system as a collection of **unravelled** atomic coordinates, that is one matter image provides us with 3xN observations.

Importantly, this means that the GPR only outputs (1+3N) terms, the tuple (E, dx1, dy1, dz1...), and that the input dimensionality of the GPR is now (1+3N).

```{code-cell} ipython3
:tags: []

ndim = Xpos_init[0].ravel().shape[0]
x_train = torch.from_numpy(
    np.row_stack([
        Xpos_init[0].ravel(),
        Xpos_init[1].ravel(),]))        
```

```{code-cell} ipython3
:tags: []

aa=np.ones(2)
bb=np.zeros(2)
np.row_stack([aa, bb])
```

```{code-cell} ipython3
:tags: []

y_train = torch.from_numpy(
    np.row_stack([
    np.append(Ye_init_np[0], Xfrcs_init[0].ravel()),
    np.append(Ye_init_np[1], Xfrcs_init[1].ravel()),
]))## should be 1, ndim+1
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
        self.base_kernel = gpytorch.kernels.RBFKernelGrad(ard_num_dims=1803) ## Important, additional dimensions
        self.covar_module = gpytorch.kernels.ScaleKernel(self.base_kernel)

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultitaskMultivariateNormal(mean_x, covar_x)

likelihood = gpytorch.likelihoods.MultitaskGaussianLikelihood(num_tasks=1804)  # Value + x-derivative + y-derivative + z-derivative
```

```{code-cell} ipython3
:tags: []

model = GPModelWithDerivativesTwo(x_train,
                               y_train,
                               likelihood)
```

```{code-cell} ipython3
:tags: []

y_train.shape
```

```{code-cell} ipython3
:tags: []

training_iter = 50
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

torch.from_numpy(Xpos_init[2].ravel()).unsqueeze(0).shape
```

```{code-cell} ipython3
:tags: []

# Set into eval mode
model.eval()
likelihood.eval()
```

```{code-cell} ipython3
:tags: []

with torch.no_grad(), gpytorch.settings.fast_computations(log_prob=False, covar_root_decomposition=False):
    test_x = torch.from_numpy(Xpos_init[2].ravel()).unsqueeze(0)
    predictions = likelihood(model(test_x))
    mean = predictions.mean
```

```{code-cell} ipython3
:tags: []

frcs = mean.squeeze(0)[1::1].reshape(-1, 3)
energy = mean.squeeze(0)[0]
```

```{code-cell} ipython3
:tags: []

mean.shape
```

```{code-cell} ipython3
:tags: []

## Average force error per atom
np.average(np.average(np.abs(frcs) - np.abs(Xfrcs_init[2]), axis=1))
```

```{code-cell} ipython3
:tags: []

np.abs(energy) - np.abs(Ye_init_np[2])
```

```{code-cell} ipython3
:tags: []

energy
```

```{code-cell} ipython3
:tags: []

predictions.mode
```

```{code-cell} ipython3
:tags: []

product.fr
```

```{code-cell} ipython3

```
