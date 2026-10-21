---
myst:
  html_meta:
    "description": "Guide to Hessian matrix calculation in eOn for vibrational analysis and hTST rate prefactors."
    "keywords": "eOn Hessian, vibrational analysis, prefactor, harmonic transition state theory, hTST"
---

# Hessian

The Hessian job is part of eOn
{cite:t}`hess-goswamiEfficientExplorationChemical2025`.
The older description of that client is Chill et al.
{cite:t}`hess-chillEONSoftwareLong2014`.

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

## How it works

eOn computes the Hessian numerically with a selectable finite-difference scheme
on the **mobile** (displaced) atoms: non-fixed atoms in `pos.con`, optionally
restricted by `[Hessian] phva_atoms` (comma-separated indices, or `All`). That
list is the hybrid/PHVA-class *active set*: atoms moved in FD. The remaining
atoms stay frozen (Li & Jensen, *Theor. Chem. Acc.* **107**, 211, 2002).

Step size is `Main.finite_difference` (historically also called finite-difference
displacement; default \(0.01\,\text{Å}\)).

### FD schemes (`[Hessian] fd_scheme`)

**one_sided** (default) — forward difference, \(\sim M\) force evaluations for
\(M = 3 N_\text{mobile}\) active DOF (plus one base gradient):

$$H_{ij} \approx -\frac{F_j(x + h e_i) - F_j(x)}{h}$$

**central** — central difference, \(\sim 2M\) force evaluations:

$$H_{ij} \approx -\frac{F_j(x + h e_i) - F_j(x - h e_i)}{2h}$$

**fourth** — fourth-order central difference, \(\sim 4M\) force evaluations
(aliases `fourth_order` and `central4`):

$$H_{ij} \approx -\frac{-F_j(x+2h e_i)+8F_j(x+h e_i)-8F_j(x-h e_i)+F_j(x-2h e_i)}{12h}$$

One-sided is preferred for classical EAM AKMC cost; central reduces \(O(h)\) bias
when validating prefactors. Fourth-order cancels the \(O(h^2)\) term as well.
The same `fd_scheme = fourth` value selects that stencil for Lanczos and
Davidson Hessian-vector products. `one_sided` and `central` leave those
products on the forward difference. The rgsaddle min-mode
(`min_mode_method = xtsci`) takes forward differences for `one_sided` and
central differences for `central` and `fourth`; see {doc}`optimizer`. Dimer
rotation, conjugate gradients, and LBFGS keep their own real steps. Complex-step is not a scheme: cutoffs and
external potentials are not holomorphic. The assembled matrix is symmetrized
\((H+H^T)/2\) before diagonalization. Mass-weighting uses
\(\tilde H_{ij} = H_{ij}/\sqrt{m_i m_j}\).

### Hessian-range coloring

Coloring is exact only when every Hessian block between two atoms farther
apart than some distance is zero. A potential reports that distance as
`Potential::hessianRange`. The value is in the same length unit as the
positions. `0`, or any non-finite value, means the range is unknown and
coloring stays off. One column is then built per coordinate.

An interaction cutoff is not this distance. A pair potential is zero past
its cutoff, so the range equals the cutoff. An atom-centred many-body
potential (EAM, Gupta, Stillinger-Weber, Tersoff) is zero only past twice
the cutoff, because the density or the bond order at an atom depends on
every neighbour inside the cutoff. A flexible water model does not have a
hard intramolecular extent, so a cutoff on the oxygen-oxygen interaction
does not make the blocks zero past any fixed distance. A message-passing
model is zero past `n_layers` times the edge cutoff, and only when that
layer count is part of the kernel.

Coloring is opt-in. `RgpotAdapter` reads `HessianRangeOf<Kernel>` and does
not look at `config().cutoff`. A kernel with a cutoff field and no
specialization stays at 0. The specialization returns the cutoff for a
pair kernel, twice the cutoff for an atom-centred many-body kernel, the
cutoff plus twice a hard intramolecular extent for a site model, or
`n_layers` times the edge cutoff for message passing. Fe-He is the
exception: its range is the larger of its pair cutoffs and twice its
largest density cutoff. Anything that is not exactly derivable returns 0.
A potential that is not an rgpot kernel opts in the same way, by
overriding `hessianRange`. An override of the deprecated
`finiteCutoff()` still supplies that distance. Everything else stays at 0.

The table of kernels eOn can wrap, at the current rgpot pin and at rgpot
`v3.5.0`, is in the developer page
{doc}`/devdocs/hessian`.

Mobile atoms are greedy-colored so that two atoms share a color only when
their closed neighborhoods of that range are disjoint. One displacement of
a whole color along x, y, or z fills every column of that color. Central
differences use two force evaluations per Cartesian direction per color,
and fourth-order differences use four, independent of the number of atoms.
Catch2 `Colored FD Hessian matches serial central difference` uses four
atoms on a line, spaced 1.2 Å, with a Lennard-Jones cutoff of 1.5 Å.
Coloring still reduces the force calls: three colors, so the central
scheme takes `1 + 2 × 3 × 3 = 19` force calls. The per-column loop takes
`1 + 2 × 12 = 25`. The colored matrix matches the serial central matrix
to `1e-8`. Catch2 `Embedding-chain coloring matches serial only at twice
the cutoff` is a planar zig-zag of eight atoms, bond 0.75 Å, with the
embedding cutoff at 1 Å. The central scheme takes 49 force calls total
and 43 once colored at `2 r_c`. The fourth-order scheme takes 97
and 85. Coloring at `r_c` takes 19 and 37. The colored matrix matches
the per-column matrix to `1e-8`. The short-range matrix misses by
0.1875 against max `|H|` 1.041667, and an off-axis entry reaches
0.541. The same central numbers hold with `remove_net_force` on and
one atom fixed, which is the path that still colors.
Catch2 `Opted-in kernels match the serial matrix at the declared range`
runs each opted-in kernel on its own zig-zag. The colored matrix matches
the per-column matrix exactly in these runs. A block past the interaction
cutoff is nonzero for every many-body kernel, and every block past the
declared range is 0. ZBL is a pair kernel, so that past-cutoff window is
empty. Its nearest-neighbour block at 1.5 Å is nonzero and every block
past `cut_global` is 0. The per-kernel call counts and the short-range
errors are in the developer page {doc}`/devdocs/hessian`.
Catch2 `Coloring falls back when every mobile atom needs its own color`
is four atoms inside a 20 Å Lennard-Jones cutoff, and the same cluster
with a batching spring of range 30 Å. Both need four colors. Both take
the column path, 25 force calls, which is one base gradient and two
central samples on 12 coordinates. The spring records a batch call, so
the colored loop did not run.
Catch2 `Fourth-order stencil matches a complex-step
oracle` evaluates `z^4 + 0.3 z^2` at `0.8`. The fourth-order value stays
within `1e-10` of the complex-step derivative. The central value is more
than `1e-3` away, and farther from the oracle than the fourth-order
value. The same matrix is built one column at a time when the range is
unknown, when `[Main] remove_net_force` couples every atom, or when a
column checkpoint is in use.

When `checkpoint_path` is empty and the potential evaluates batches, the
displaced structures go out together, at most 32 columns per call.
`RgpotPot` takes that path once `[RgpotPot] ranks_per_image` has more than
one calculator group. Structure `j` runs on group `j` modulo the group
count, and every rank receives every result. A prefactor Hessian uses the
same call.

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
# fd_scheme = fourth
# resume = true
# checkpoint_path = hessian.ckpt
zero_freq_value = 1e-6
```

The output `results.dat` records force-call counts; `hessian.dat` contains the
mass-weighted matrix. `[Main] quiet = true` suppresses the write log line only.
Eigenvalues (squared frequencies) are obtained by diagonalizing the symmetrized
matrix (ColMajor eigen solve in the client).

The job drops each eigenvalue at or below `zero_freq_value` (default
`1e-6`). A free cluster has 6 of those eigenvalues, a linear cluster has
5, a periodic cell has 3, and a structure with any fully fixed atom has
none. The job still drops a different count, and it logs a warning. It
does not stop for that warning.

With `[Hessian] write_modes = true` (the default), `modes.con` holds one frame
per eigenvalue, lowest first. Each frame is the input structure with the mode
in its `displacements` section: the Cartesian displacement \(x = q/\sqrt{m}\)
of the mass-weighted eigenvector \(q\), scaled to unit norm, and zero on atoms
outside the mobile set. The frame metadata carries `mode_eigenvalue` in
eV/(Å² amu), `hbar_omega` in eV and `wavenumber` in cm⁻¹. The last two are
negative for an imaginary mode. Any readcon reader gets the mode as an
N x 3 array beside the positions:

```{code-block} python
import readcon

for frame in readcon.read_con("modes.con"):
    print(frame.metadata["wavenumber"], frame.disp.shape)
```

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

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: HS_
keyprefix: hess-
---
```
