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

eOn computes the Hessian numerically with a selectable finite-difference scheme
on the **mobile** (displaced) atoms: non-fixed atoms in `pos.con`, optionally
restricted by `[Hessian] atom_list` (comma-separated indices, or `All`). That
list is the hybrid/PHVA-class *active set* — atoms that are moved in FD, not
the frozen environment (Li & Jensen, *Theor. Chem. Acc.* **107**, 211, 2002).

Step size is `Main.finite_difference` (historically also called finite-difference
displacement; default \(0.01\,\text{Å}\)).

### FD schemes (`[Hessian] fd_scheme`)

**one_sided** (default) — forward difference, \(\sim M\) force evaluations for
\(M = 3 N_\text{mobile}\) active DOF (plus one base gradient):

$$H_{ij} \approx -\frac{F_j(x + h e_i) - F_j(x)}{h}$$

**central** — central difference, \(\sim 2M\) force evaluations:

$$H_{ij} \approx -\frac{F_j(x + h e_i) - F_j(x - h e_i)}{2h}$$

One-sided is preferred for classical EAM AKMC cost; central reduces \(O(h)\) bias
when validating prefactors. The assembled matrix is symmetrized \((H+H^T)/2\)
before diagonalization. Mass-weighting uses \(\tilde H_{ij} = H_{ij}/\sqrt{m_i m_j}\).

### Column resume

Long partial Hessians can checkpoint FD columns to `checkpoint_path` (e.g.
`hessian.ckpt`) with `resume = true`. Interrupted jobs continue from the next
column; a successful run removes the checkpoint.

## Usage

```{code-block} ini
[Main]
job = hessian
finite_difference = 0.01

[Hessian]
atom_list = All
fd_scheme = one_sided
# fd_scheme = central
# resume = true
# checkpoint_path = hessian.ckpt
zero_freq_value = 1e-6
```

The output `results.dat` records force-call counts; `hessian.dat` holds the
mass-weighted matrix when `quiet = false`. Eigenvalues (squared frequencies)
are obtained by diagonalizing the symmetrized matrix (ColMajor eigen solve in
the client).
