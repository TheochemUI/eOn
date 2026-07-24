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
restricted by `[Hessian] phva_atoms` (comma-separated indices, or `All`). That
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
phva_atoms = All
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

## Free/fixed versus active (PHVA)

**free/fixed** on `pos.con` / `Matter` is the optimizer mask: which atoms may
move during geometry steps. Under an active-volume setup that set can be large
(the entire movable shell).

**`phva_atoms`** is the PHVA *active* set (Li & Jensen, *Theor. Chem. Acc.*
**107**, 211, 2002): atoms that are *displaced* when building a partial
Hessian or a matrix-free Hessian-vector product. It is always intersected with
free flags so frozen atoms are never moved. Default `All` means
active = free (historical behavior).

The same polarity applies to matrix-free min-mode. Free/fixed is left unchanged;
an explicit list sets the Krylov / Ritz dimension to \(3 N_\mathrm{active}\)
rather than \(3 N_\mathrm{free}\):

```{code-block} ini
[Lanczos]
phva_atoms = 12, 15, 18, 21
tolerance = 0.01
max_iterations = 20

[Davidson]
phva_atoms = 12, 15, 18, 21
```

In pyeonclient:

```{code-block} python
import pyeonclient as pyec

mobile = pyec.resolve_mobile_atoms(matter, "12,15,18,21")
lanczos.compute(matter, direction, atoms=mobile)
# or: params.lanczos_phva_atoms = "12,15,18,21"
#     lanczos.compute(matter, direction)
```
