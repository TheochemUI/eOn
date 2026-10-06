---
myst:
  html_meta:
    "description": "eOn rgpot potential: in-process rgpot NWChem and CPMD, including a CPMDParams file and MPI calculator groups."
    "keywords": "eOn, rgpot, CPMD, libcpmdc, params_path, ranks_per_image"
---

# RgpotPot (direct in-process rgpot)

```{versionadded} 2.16.0
```

On non-Windows builds, potential `rgpot` links
[rgpot](https://github.com/OmniPotentRPC/rgpot) and loads `libnwchemc.so` or
`libcpmdc.so` with `dlopen` in the eOn process.
The ini name is `rgpot`.
The cpmdc engine is documented at [cpmdc.rgoswami.me](https://cpmdc.rgoswami.me).

Potserv clients, `eonclient --serve`
([Serve mode](project:serve_mode.md)), and SocketNWChem are separate.
Overview: [rgpot integration](project:rgpot_integration.md).

## Build

```{code-block} bash
meson setup bbdir -Dwith_tests=true
meson compile -C bbdir
```

The build needs Cap'n Proto headers and libraries. Method parameters are
Cap'n Proto messages passed into the C ABI. The `rgpot` Meson subproject is
`subprojects/rgpot.wrap`, and the build pulls `nwchempot_dep` and
`cpmdpot_dep`. Serve mode uses `ptlrpc_dep` under `-Dwith_serve`.

`dependency('rgpot')` prefers an installed rgpot. `cpmdc_bind_calculator`
and `rgpot::bindCalculators` first ship in rgpot 3.3.0. This build requires
3.4.0, which adds `rgpot::calculatorsUseMpi`. An installed 3.3.0 fails that
version check. The Meson wrap is the fallback. A CPMD run that splits ranks
uses the flag in the CPMD section.

## Configuration

`eon/config.yaml` and the model list the same `[RgpotPot]` keys. `eonclient`
reads those keys from `config.ini`.

```{code-block} ini
[RgpotPot]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.RgpotPot
```

### NWChem

```{code-block} ini
[Main]
job = point

[Potential]
potential = rgpot

[RgpotPot]
backend = nwchemc
basis = sto-3g
theory = scf
scf_type = rhf
charge = 0
multiplicity = 1
```

`engine_path`, or `engine_library` when `engine_path` is empty, selects
`libnwchemc.so`. For this backend the client also reads `NWCHEMC_LIBRARY`
and `RGPOT_NWCHEMC_ENGINE`.

`input_block`, or `RGPOT_NWCHEM_INPUT_BLOCK` when the key is empty, is
NWChem `inputBlocks` text. When `theory` is `dft` and `scf_type` names an
exchange-correlation functional such as `b3lyp`, the client writes a short
density functional theory (DFT) block.

### Metatomic

```{code-block} ini
[Potential]
potential = rgpot

[RgpotPot]
backend = metatomic
model_path = /path/to/model.pt
device = cpu
```

`engine_path` may point at `libmetatomic_engine.so`. The client also reads
`RGPOT_METATOMIC_ENGINE` and `METATOMIC_ENGINE`.

### xTB

This backend leaves `-Dwith_xtb=false` on the eOn build. The engine is
`libxtb_engine.so`, loaded at run time.

```{code-block} ini
[Potential]
potential = rgpot

[RgpotPot]
backend = xtb
paramset = GFN2xTB
accuracy = 1.0
```

`engine_path` may point at `libxtb_engine.so`. The client also reads
`RGPOT_XTB_ENGINE` and `XTB_ENGINE`.

### UMA

The engine is `libuma_engine.so`, loaded at run time through the generic
engine interface. `model_path` names the ahead-of-time compiled model,
`task_name` selects its task head (default `omol`), and `charge` and
`multiplicity` set the total charge and spin.

```{code-block} ini
[Potential]
potential = rgpot

[RgpotPot]
backend = uma
model_path = /path/to/uma.pt2
task_name = omol
device = cpu
```

`engine_path` may point at `libuma_engine.so`. The client also reads
`RGPOT_UMA_ENGINE`. A molecule whose `.con` file carries no cell gets a
25 angstrom diagonal box.

## CPMD

`backend = cpmdc` runs one Car-Parrinello molecular dynamics (CPMD)
session inside `eonclient`. The shared library is
[libcpmdc](https://github.com/OmniPotentRPC/cpmdc).
The engine documentation is at [cpmdc.rgoswami.me](https://cpmdc.rgoswami.me).

Calculator groups need rgpot built with the Message Passing Interface (MPI).
Pass the flag to the wrap:

```{code-block} bash
meson setup bbdir -Drgpot:with_mpi=enabled
meson compile -C bbdir
```

The flag configures the wrap. A `pkg-config` rgpot must already be an MPI
build. If it is not, `ranks_per_image` greater than 0 raises.

Neither `librgpot` nor eOn's `librgpot_pot` links MPI, so an `eonclient`
started without `mpirun` loads no MPI library and pays no MPI start-up
cost. The MPI side is loaded on demand: rgpot's `librgpot_mpi`, and
eOn's `librgpot_pot_mpi`, which holds the construction agreement and the
grouped exit handler below. eOn loads it only under an MPI launcher,
from `EON_RGPOT_MPI_LIBRARY`, the directory of `librgpot_pot`, or the
linker path. It builds when rgpot has MPI: the wrap with
`-Drgpot:with_mpi=enabled`, or an installed rgpot whose `rgpot.pc`
defines `RGPOT_HAS_MPI`; eOn's `-Dwith_mpi=enabled` also builds it.
After a failed engine call rgpot asks for `MPI_Abort` at exit, and the
exit handler aborts the world instead of waiting in `MPI_Finalize` for
ranks left inside a CPMD collective.

eOn's `-Dwith_mpi=enabled` option builds the client/server program. Calculator
groups are this page's launch, `mpirun -np N eonclient`, with rgpot built
`-Drgpot:with_mpi=enabled`.

### CPMDParams file and the cpmd section

Three layers build the message CPMD receives.

`params_path` on `[RgpotPot]` is a Cap'n Proto CPMDParams file. It owns
`functional`, `cutOffRy`, `charge`, `multiplicity`, `title`, `memoryMb`,
`inputSections`, and `inputBlocks`. `RGPOT_PARAMS_PATH` overrides the
ini key. Scalar keys are not written over a file that loaded. Every rank
reads the file before `MPI_Comm_split`. When any rank cannot read it,
every rank throws that error and the split does not run.

When `params_path` is empty, `[cpmd]` supplies those scalars and overrides
the copies on `[RgpotPot]`. Three spellings set the cutoff on either
section: `cutOffRy`, then `cutoff_ry`, then `cpmd_cut_off_ry`. The first
one present wins. `functional` wins over `cpmd_functional`. The client and
the Python server accept every spelling; the `Cpmd` and `RgpotPot` models
write only `cutOffRy` on `[cpmd]` and `cutoff_ry` on `[RgpotPot]`.
`[cpmd]` is read only when `backend` is `cpmd`, `cpmdc`, or `cpmdpot`.

`engine_path`, `engine_library`, `engine_root`, `scratch_dir`,
`permanent_dir`, and `ranks_per_image` stay on `[RgpotPot]`. They place
the process and apply after either source.

`input_block` is one more `inputBlocks` entry. `[cpmd]` supplies it, then
`[RgpotPot]`, then `RGPOT_CPMD_INPUT_BLOCK` when those strings are empty.
Sections loaded from the file stay. A block already in the file stays,
and the ini text follows it. An empty `input_block` leaves the file's
blocks unchanged.

```{code-block} ini
[Potential]
potential = rgpot

[RgpotPot]
backend = cpmdc

[cpmd]
functional = BLYP
cutOffRy = 70.0
charge = 0
multiplicity = 1
```

```{eval-rst}
.. autopydantic_model:: eon.schema.Cpmd
```

That example is the scalar layer. The point, minimization, and band
examples below keep `params_path`, which is the file layer.

Write the message as Cap'n Proto text. The field names are in the
[write-cpmdparams how-to](https://github.com/OmniPotentRPC/cpmdc/blob/main/docs/source/howto/write-cpmdparams.rst).
From a cpmdc checkout, encode it:

```{code-block} bash
capnp encode schema/Potentials.capnp CPMDParams \
  < cluster.params.txt > cluster.params.bin
```

`capnp encode` writes the flat message `eonclient` reads. A field name
that is not in the schema fails at this step.

### Environment

```{code-block} bash
export CPMDC_LIBRARY=/path/to/libcpmdc.so
export CPMDC_PSEUDO_DIR=/path/to/PP_LIBRARY
```

`RGPOT_CPMDC_ENGINE` is the other library path the client reads. rgpot
then tries `RGPOT_CPMD_ENGINE`, then `libcpmdc.so` on the loader path.

libcpmdc reads `CPMDC_PSEUDO_DIR` for the pseudopotential directory. When
that variable is unset, it reads `CPMD_PP_LIBRARY_PATH`.

`RGPOT_CPMD_INPUT_BLOCK` fills `input_block` when the `[cpmd]` key and
the `[RgpotPot]` key are empty. The text is appended to `inputBlocks`.
cpmdc places those blocks ahead of the sections it generates. Sections
in the CPMDParams file are not removed.

`permanent_dir` is the CPMD `FILEPATH` for `RESTART` files. `scratch_dir`
is the fallback directory.

### Files from libcpmdc

A library force call writes no `RESTART.1`, `LATEST`, `GEOMETRY`, or
`GEOMETRY.xyz`. The orbitals for the next call stay in memory. `cpmd.x`
writes those four files in `permanent_dir` when that key is set, otherwise
in `scratch_dir`, otherwise in the working directory.

eOn's own files are written by rank 0. The names depend on the job.

### A point

The geometry is `pos.con`. `ranks_per_image = 0` is one group of all
ranks. One process is that group.

```{code-block} ini
[Main]
job = point

[Potential]
potential = rgpot

[RgpotPot]
backend = cpmdc
params_path = cluster.params.bin
ranks_per_image = 0
```

```{code-block} bash
eonclient
```

`results.dat` records the energy in eV, the maximum force in eV/A, the
force-call count, and the termination status.

`mpirun -np 4 eonclient` with the same `ranks_per_image` is still one
group. Those 4 ranks share one CPMD session. Rank 0 runs the job. The
other ranks serve force requests.

### A minimization

The input geometry is `pos.con`. The job writes `min.con` and
`results.dat`. `min.con` is the minimized structure. `results.dat` records
the energy in eV, the force-call count, and whether the run converged.

`converged_force` is a threshold in eV/A.

```{code-block} ini
[Main]
job = minimization

[Potential]
potential = rgpot

[Optimizer]
opt_method = lbfgs
converged_force = 0.01
max_iterations = 1000

[RgpotPot]
backend = cpmdc
params_path = cluster.params.bin
ranks_per_image = 0
```

Launch it with the same `eonclient` or `mpirun` line as the point job.

### A 7-image nudged elastic band (NEB)

`images = 7` asks for 7 intermediate images. A band update evaluates
each image whose positions changed.

The reactant and the product are one batch at the start. The reactant
runs on group 0. The product runs on group 1.

Seven groups of 4 ranks need 28 ranks. Intermediate image 1 runs on
group 0, and image 7 runs on group 6. Every image keeps its group from
one iteration to the next, because the group holds that image's
orbitals.

A batch puts at most ceil(M / G) of its M systems on one group. An
update of only some images, which would stack them on the groups that
own them, moves the excess to the least-loaded groups, and a moved image
stays on its new group. Groups that own fewer systems than the busiest
group repeat their last system into scratch, and an empty group repeats
system 0 of that batch, so every group enters the engine the same number
of times. A repeat of a group's last system is the stored result of its
session and runs no SCF. The systems count omits those repeats.

With `ci_mmf = true` the climbing image takes improved-dimer steps on its
own. Each step evaluates the moved centre and its forward image as one
batch. A rotation trial depends on the previous one and runs alone:
every group evaluates that one structure, and the run keeps the result
of the group that last evaluated the geometry nearest to it. The usage
line counts the call on that group. A saddle search does the same, and
with `min_mode_method = lanczos` the second system of that batch is the
first Krylov product's displaced image.

A Hessian or a prefactor with an empty `checkpoint_path` sends its
displaced structures through these groups the same way.

Each group keeps the converged orbitals of every system it evaluates,
under that system's key: the image index for a band, the bead index for a
ring polymer, the position in the batch otherwise, and one key for every
single request. The next SCF of an image or a bead starts from its own
orbitals of the previous step rather than from those of whichever system
the group evaluated last. One group (`ranks_per_image = 0`) batches as
well, so the keys hold there too. The engine needs
`cpmdc_session_select_orbitals`; an older libcpmdc keeps one stored copy
per group, and so does an rgpot without `CPMDPot::selectOrbitals` (3.4.0
and older).

### Choosing the number of groups

A batch of M systems on G groups takes ceil(M / G) rounds, so its load
balance is M / (G ceil(M / G)). Pick G to divide the systems of the
batch that repeats: the interior images of a band (`images`), or the
beads of a ring (`path_beads`, `pi_beads`). Seven images run evenly on 1
or 7 groups; on 2 groups the second round leaves one group idle and the
load balance is 7/8. Sixteen beads run evenly on 1, 2, 4, 8 or 16 groups.
Past that, fewer groups of more ranks often win: a CPMD SCF step of a
small cell scales over the ranks of one group better than the groups
share a node. For the Si3N4 Geo1 cell of this page, one group of 28
ranks took 141.7 s for 9 force calls against 160.6 s for 7 groups of 4.
Those figures are wall times. The minimized energy of that cell is
`potential_energy` in `results.dat`, and this page does not quote it.
Keep to one rank per physical core.

The run reads `reactant.con` and `product.con`.

```{code-block} ini
[Main]
job = nudged_elastic_band

[Potential]
potential = rgpot

[Optimizer]
opt_method = lbfgs
converged_force = 0.01
max_iterations = 1000

[Nudged Elastic Band]
images = 7
converged_force = 0.01

[RgpotPot]
backend = cpmdc
params_path = cluster.params.bin
ranks_per_image = 4
```

```{code-block} bash
mpirun -np 28 eonclient
```

Rank 0 runs the job and writes the files. The other ranks serve force
requests. When rank 0 stops them it prints how the job used the groups:
the systems each group evaluated, its seconds inside engine calls, and
its idle share of the driver's wall time in grouped requests. The last
line gives the POP ratios. Load balance is the mean busy time over the
largest. Communication efficiency is the largest busy time over that
wall time, so a wait for the slowest image of a batch, a serial request
on group 0, or a result broadcast lowers it. Parallel efficiency is the
product of the two.

For the Si3N4 band of this page started from a nearly converged path,
on 8 ranks as 2 groups of 4:

```{code-block} text
RgpotPot: calculator groups: 2 groups, 10 batches and 6 single requests, 1216.8 s grouped wall (109.3 s single)
 group  systems     busy_s   idle
     0       34    1187.44   2.4%
     1       22     889.36  26.9%
load balance 0.874, communication efficiency 0.976, parallel efficiency 0.853
```

At exit, rank 0 sends a stop. Every rank calls `cpmdc_finalize`
while MPI is still up, then `MPI_Finalize`, then `_Exit`. `_Exit` returns
to the kernel, and the dynamic linker does not run the CPMD or MPI
library destructors on a finalized world. Worker ranks use status 0.
Rank 0 uses the job status, and a calculator error on another rank is
part of the message it prints.

The job writes `results.dat`, `neb.con`, `neb.dat`, and `sp.con`.
`sp.con` is the highest-energy image. A spline maximum more than 0.05 eV
above the reactant also produces `peak00_pos.con` and `peak00_mode.dat`,
numbered from `00`. The default of `setup_mmf_peaks` leaves those files on.

## SocketNWChem

SocketNWChem speaks the i-PI socket, and eOn listens. The direct potential
(`rgpot`) loads `libnwchemc.so` or `libcpmdc.so` in the eOn process.
SocketNWChem keeps an external NWChem process warm across calls. `rgpot`
keeps one session, and that session keeps the orbitals.

## Implementation notes

`RgpotPot` owns one translation unit. That unit includes rgpot headers and
no eOn potential header. The Cap'n Proto name `Potential` stays off eOn's
`Potential` class.

Energies are in eV and forces are in eV/A. rgpot converts the Hartree and
Hartree/bohr values that the engine returns.
