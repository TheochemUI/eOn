---
myst:
  html_meta:
    "keywords": "eOn IRC, intrinsic reaction coordinate, saddle search, rgsaddle, Morse Pt, connected minima"
    "description": "Check which minima a converged eOn saddle connects by following its intrinsic reaction coordinate with rgsaddle."
---

# Follow the intrinsic reaction coordinate from a saddle

This how-to shows you how to check which two minima a saddle connects.
It follows the intrinsic reaction coordinate (IRC) {cite:p}`irc-fukuiPathChemicalReactions1981` downhill from the saddle in both directions.
You add three lines to `[Saddle Search]` and run `eonclient` once.
You get the path, its energy profile and the two relaxed end states.

The example is the Morse Pt saddle from {doc}`pt_saddle`: one free Pt atom on a fixed Pt slab.
The log line, the `results.dat` lines and the figure come from one run.
Its files sit in `docs/source/tutorials/systems/data/irc_pt/`.
Other runs give the table of step sizes and the cost section, and each says which runs.

## Before you start

You need an `eonclient` built with rgsaddle.
From the root of the eOn checkout:

```{code-block} bash
meson setup bbdir --buildtype release -Dwith_rgsaddle=true
meson compile -C bbdir
```

Meson finds rgsaddle through pkg-config.
If you installed rgsaddle under your own prefix, add that prefix's `lib/pkgconfig` to `PKG_CONFIG_PATH` first.
The continuous integration (CI) job in `.github/workflows/ci_rg_backends.yml` builds and installs the pinned rgsaddle this way.

The IRC needs rgsaddle C application binary interface (ABI) minor 14 or newer.
Both the build and the library loaded at run time must meet it.
eOn refuses `irc_after = true` before the saddle search starts when either is older.
A missing feature never costs you a saddle search.

## Set up the job

Copy the inputs into a new directory.
Run these commands from the root of the eOn checkout:

```{code-block} bash
mkdir irc_pt
cp docs/source/tutorials/systems/data/pt_saddle/pos.con irc_pt/
cp docs/source/tutorials/systems/data/pt_saddle/displacement.con irc_pt/
cp docs/source/tutorials/systems/data/irc_pt/direction.dat irc_pt/
cp docs/source/tutorials/systems/data/irc_pt/config.ini irc_pt/
cd irc_pt
```

`direction.dat` gives the dimer its starting direction.
Without it this search ends on the reactant with a zero eigenvalue.
eOn then skips the IRC with a warning.
A point with no negative eigenvalue is not a saddle.

This is the whole `config.ini`:

```{literalinclude} data/irc_pt/config.ini
:language: ini
```

The last three lines are the only change from a plain saddle search.
`irc_after = true` turns the IRC on.
`irc_step_size` is the arc length of one step in mass-weighted coordinates, in amu{sup}`1/2` Angstroms.
Here amu is the atomic mass unit.
For a Pt atom (195 amu) a step of 0.5 moves the atom about 0.036 Angstroms.
`irc_converged_force` is the largest force on any atom at which a branch counts as having reached a minimum.
It is in electronvolts/Angstrom (eV/Angstrom).
This page sets it to 0.01, below the default of 0.05.
The section on choosing the step size and force threshold below shows why.

## Run it

Run the client you built, from inside `irc_pt`:

```{code-block} bash
../bbdir/client/eonclient
```

The log ends with one line for the IRC:

```{code-block} text
[Saddle Search] IRC: 2 of 2 branches at a minimum, 228 force calls, 2 more to relax the ends
```

`results.dat` gains seven lines after the usual saddle entries.
These are the lines from the run, and the whole file is {download}`data/irc_pt/results.dat`:

```{code-block} text
-1462.008706 potential_energy_saddle
-1.010010 final_eigenvalue
-1462.166782 potential_energy_reactant
2 irc_branches_at_minimum
1 irc_ends_relaxed
1 irc_connects_reactant
228 force_calls_irc
2 force_calls_irc_minimization
-1462.166783 potential_energy_irc_min1
-1462.179029 potential_energy_irc_min2
```

Read them in this order:

1. `irc_connects_reactant = 1` says one relaxed end matched `pos.con` under the `[Structure Comparison]` settings.
   This is the check that the saddle belongs to the state you started from.
2. `potential_energy_irc_min1` and `potential_energy_irc_min2` are the two relaxed ends.
   Here min1 is the reactant and min2 is a new minimum 0.012 eV lower.
3. `irc_branches_at_minimum` counts the branches that ended at a minimum within `irc_max_steps`.
   eOn counts a branch when rgsaddle reports a minimum.
   rgsaddle reports one only after the first step and where the largest atomic force is at or below `irc_converged_force`.
   A branch that stops short still gets relaxed.
   A count of 2 does not by itself mean the path is complete.
   The first row of the step-size table is a run at the defaults.
   There both branches count, but the min2 branch ends 0.0152 Angstroms from the saddle.
4. `total_force_calls` is 270 here.
   It counts 1 call for the energy of `pos.con`, 39 for the saddle, 228 for the path and 2 for the ends.
   A run of the same job with the option off takes 40.

The job also writes four files:

| File | Contents |
|------|----------|
| `irc.con` | the path as frames, from the min1 end through the saddle to the min2 end, with each frame's energy in its header |
| `irc.dat` | one row per frame, laid out like `neb.dat` with `max_force` as the fourth column: `img`, `rxn_coord`, `energy`, `max_force` |
| `irc_min1.con` | the relaxed min1 end, against the mode |
| `irc_min2.con` | the relaxed min2 end, along the mode |

In `irc.dat`, `rxn_coord` is the Cartesian distance along the path in Angstroms.
It is negative on the min1 side and zero at the saddle.
`energy` is relative to the saddle.
min1 and min2 use the same sides as a process search, which relaxes from the saddle minus and plus the mode.

## Plot the energy along the path

`plot_irc.py` in the data directory plots `irc.dat` and marks the relaxed ends from `results.dat`.
It needs NumPy and Matplotlib.
From inside `irc_pt`:

```{code-block} bash
python ../docs/source/tutorials/systems/data/irc_pt/plot_irc.py . irc_energy.png
```

```{figure} data/irc_pt/irc_energy.png
:alt: Energy relative to the saddle against distance along the IRC, where both branches fall from the saddle at zero, the min1 branch ends near -0.158 eV at -0.79 Angstroms and the min2 branch near -0.170 eV at 0.81 Angstroms.

Energy along the IRC from the Morse Pt saddle.
Open circles are the ends after relaxation.
The min1 branch is densely sampled
because rgsaddle halves the step after each rejected trial
and keeps the shorter step for the rest of that branch.
The plot reads {download}`data/irc_pt/irc.dat` and {download}`data/irc_pt/results.dat` from this run.
```

## Choose the step size and force threshold

The defaults are `irc_step_size = 0.1` and `irc_converged_force = 0.05`.
The threshold is the path stop of Sella's IRC and the default of the rgsaddle session.
The table lists six runs of this job from the same saddle, with the rgsaddle that the CI job pins.
They differ only in these two keys, and the last row is the run in `data/irc_pt/`.
The end of a branch is its last IRC point, before eOn relaxes it.
The end columns give its `rxn_coord` in `irc.dat`.

| `irc_step_size` | `irc_converged_force` | path force calls | end force calls | min1 end (Angstroms) | min2 end (Angstroms) |
|---|---|---|---|---|---|
| 0.1 (default) | 0.05 (default) | 123 | 22 | -0.788 | 0.0152 |
| 0.1 (default) | 0.01 | 1021 | 2 | -0.788 | 0.814 |
| 0.2 | 0.05 (default) | 234 | 2 | -0.788 | 0.814 |
| 0.2 | 0.01 | 234 | 2 | -0.788 | 0.814 |
| 0.5 | 0.05 (default) | 58 | 23 | -0.0402 | 0.814 |
| 0.5 | 0.01 | 228 | 2 | -0.788 | 0.814 |

All six runs relax into the same two minima.
At the defaults the min2 branch stops 0.0152 Angstroms from the saddle, where the largest force is 0.017 eV/Angstrom.
At step 0.5 and a threshold of 0.05 the min1 branch stops 0.0402 Angstroms out, where it is 0.049 eV/Angstrom.
Both forces are below 0.05 eV/Angstrom.
The pinned rgsaddle reports these points as minima, so eOn counts those branches.
A step of 0.2 traces both branches at either threshold.
At a threshold of 0.01 every step traces both branches.
Of the runs that trace both, the one at 0.5 and 0.01 takes the fewest force calls, so this page uses it.

If you only need the two minima, the default threshold is cheaper here.
The path and the ends take 145 force calls at the default step and 81 at step 0.5, against 230 for the run on this page.

`irc_max_steps` (default 1000) caps the steps per branch.

## What it costs

A process search already checks a saddle with two minimizations displaced along the mode.
The IRC does not replace that cheaper check.
We tried both on four systems.
These were the Morse Pt adatom saddle of this page and the 7- and 13-atom Lennard-Jones clusters.
The fourth was the Pt island of `examples/akmc-pt`.
The first three ran at steps from 0.05 to 0.5 and at thresholds of 0.05 and 0.01, and the island at step 0.1.
Both reached the same two minima every time, and the IRC always used more force calls.
Use the IRC when you need the path.
Also use it when you suspect that the displaced minimizations slide into a minimum the saddle does not connect.

rgsaddle keeps a dense quasi-Newton model of the free coordinates and diagonalizes it more than once per step.
The Pt island has 175 free atoms.
At step 0.1 its 14 IRC steps took about 120 s for 16 force calls, at either threshold.
The two runs took 120.00 s at 0.05 and 118.15 s at 0.01.
Both branches stopped about one step from the saddle.
Relaxing the two ends took 0.08 s for 118 force calls.
Keep the IRC to small free regions for now.

## See also

- Every `[Saddle Search]` key: {doc}`/user_guide/saddle_search`
- The saddle search on its own: {doc}`pt_saddle`
- The displaced minimizations of a process search: {doc}`/user_guide/process_search`

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: IRC_
keyprefix: irc-
---
```
