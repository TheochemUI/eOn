---
myst:
  html_meta:
    "description": "Guide to Hessian matrix calculation in eOn for vibrational analysis and hTST rate prefactors."
    "keywords": "eOn Hessian, vibrational analysis, prefactor, harmonic transition state theory, hTST"
---

# Hessian

The Hessian matrix (second derivative of the potential energy with respect to
atomic coordinates) is used for:

- **Vibrational frequency analysis**: Eigenvalues of the mass-weighted Hessian
  give squared vibrational frequencies. Positive eigenvalues correspond to stable
  modes; negative eigenvalues indicate saddle point character.
- **hTST prefactors**: The harmonic transition state theory rate constant
  requires the product of frequencies at the minimum and saddle point
  (see [prefactor](project:prefactor.md)).
- **Saddle verification**: A first-order saddle point has exactly one negative
  Hessian eigenvalue.

## How It Works

eOn computes the Hessian numerically via central finite differences:

$$H_{ij} = \frac{F_i(x + \delta e_j) - F_i(x - \delta e_j)}{2\delta}$$

where $\delta$ is the finite difference step size (`min_displacement`). This
requires $2 \times 3N$ gradient evaluations for $N$ atoms (or $2 \times 3N_\text{free}$
if some atoms are frozen).

The Hessian is then mass-weighted:

$$\tilde{H}_{ij} = \frac{H_{ij}}{\sqrt{m_i m_j}}$$

and diagonalized to obtain vibrational frequencies $\nu_k = \sqrt{\lambda_k} / (2\pi)$.

## Usage

The Hessian job computes and reports the vibrational frequencies:

```{code-block} ini
[Main]
job = hessian

[Hessian]
min_displacement = 0.001
```

The output `results.dat` contains the eigenvalues (squared frequencies) of the
mass-weighted Hessian matrix.

## Configuration

```{code-block} ini
[Hessian]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.HessianConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: HESS_
keyprefix: hess-
---
```
