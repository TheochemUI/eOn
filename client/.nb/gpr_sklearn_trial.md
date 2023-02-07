---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.4
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

```python tags=[]
import numpy as np
import numpy.typing as npt
import sklearn

import matplotlib.pyplot as plt
import copy

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, Matern

import pyeonclient as ec
```

```python tags=[]
params = ec.Parameters()
params.load("../../examples/neb-al/config.ini") # NEB-Al
ec.log_init(params, "blah.log") # Just in case
```

```python tags=[]
product = ec.Matter(params)
product.con2matter("../../examples/neb-al/product.con")
reactant = ec.Matter(params)
reactant.con2matter("../../examples/neb-al/reactant.con")
```

```python tags=[]
neb = ec.NudgedElasticBand(reactant, product, params)
nmatter = neb.neb_images
```

```python tags=[]
X = [x.positions.ravel() for x in nmatter]
```

```python tags=[]
y = [x.pot_energy for x in nmatter]
```

```python tags=[]
Yf = [x.forces.ravel() for x in nmatter]
```

```python tags=[]
np.concatenate([Yf], axis=1 ).shape
```

```python tags=[]
np.concatenate([X], axis=1 ).shape
```

```python tags=[]
#kernel = 1 * RBF(length_scale=1.0, length_scale_bounds=(1e-2, 1e2))
#gaussian_process = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=9)
noise_std = 0.5
kernel = 1 * RBF(length_scale=1.0, length_scale_bounds=(1e-2, 1e9))
#kernel = 1.0 * Matern(length_scale=[1.0]*X[0].shape[0], nu=1.5)
#gaussian_process = GaussianProcessRegressor(kernel=kernel, random_state=0)
gaussian_process = GaussianProcessRegressor(
    kernel=kernel, alpha=noise_std**2, n_restarts_optimizer=13
)
gaussian_process.fit(np.concatenate([X], axis=1 ), np.concatenate([Yf], axis=1))
gaussian_process.kernel_
```

```python tags=[]
gprp = lambda x: gaussian_process.predict(x.reshape(1, -1), return_std=True)
```

```python tags=[]
neb.compute()
```

```python tags=[]
nm_fin = neb.neb_images
```

```python tags=[]
[x.pot_energy for x in nm_fin]
```

```python tags=[]
[gprp(x.positions.ravel()) for x in nm_fin]
```

```python tags=[]
[gprp(x.ravel()) for x in np.concatenate([X], axis=1)]
```

```python tags=[]
y
```

## SKL Baseline

```python tags=[]
X = np.linspace(start=0, stop=10, num=1_000).reshape(-1, 1)
y = np.squeeze(X * np.sin(X))
```

```python tags=[]
plt.plot(X, y, label=r"$f(x) = x \sin(x)$", linestyle="dotted")
plt.legend()
plt.xlabel("$x$")
plt.ylabel("$f(x)$")
_ = plt.title("True generative process")
```

```python tags=[]
rng = np.random.RandomState(1)
training_indices = rng.choice(np.arange(y.size), size=6, replace=False)
X_train, y_train = X[training_indices], y[training_indices]
```

```python tags=[]
kernel = 1 * RBF(length_scale=1.0, length_scale_bounds=(1e-2, 1e2))
gaussian_process = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=9)
gaussian_process.fit(X_train, y_train)
gaussian_process.kernel_
```

```python tags=[]
mean_prediction, std_prediction = gaussian_process.predict(X, return_std=True)

plt.plot(X, y, label=r"$f(x) = x \sin(x)$", linestyle="dotted")
plt.scatter(X_train, y_train, label="Observations")
plt.plot(X, mean_prediction, label="Mean prediction")
plt.fill_between(
    X.ravel(),
    mean_prediction - 1.96 * std_prediction,
    mean_prediction + 1.96 * std_prediction,
    alpha=0.5,
    label=r"95% confidence interval",
)
plt.legend()
plt.xlabel("$x$")
plt.ylabel("$f(x)$")
_ = plt.title("Gaussian process regression on noise-free dataset")
```

```python

```
