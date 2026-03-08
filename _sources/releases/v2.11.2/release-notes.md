---
myst:
  html_meta:
    "description": "Detailed release notes for eOn v2.11.2, covering metatomic FPE fix, uncertainty default change, and improved dimer log alignment."
    "keywords": "eOn release notes, metatomic, FPE, SIGFPE, uncertainty, dimer, IDimerRot"
---

# Release notes

## [v2.11.2] - 2026-03-02

### Fixed

#### Metatomic SIGFPE with in-process libtorch

The built-in `metatomic` potential triggered `SIGFPE: Invalid operation` during
libtorch forward/backward passes because eOn enables global FPE trapping at
startup.  Libtorch internals (SiLU via sleef, autograd intermediates) produce
benign NaN/Inf that get trapped.  The constructor and `force()` method now
suppress FPE trapping for the duration of torch operations using `FPEHandler`,
matching the pattern already used by ASE_ORCA, ASE_NWCHEM, and AtomicGPDimer.

#### Uncertainty threshold default too eager

`uncertainty_threshold` defaulted to `0.1` in both C++ and Python.  Most models
lack uncertainty outputs, so the default triggered a noisy exception+catch in the
metatomic constructor for no benefit.  The default is now `-1` (disabled),
making uncertainty checking opt-in.

#### IDimerRot column misalignment

The `[IDimerRot]` force column placeholder was 10 dashes but `[Dimer]` uses
`{:18.5e}` (18 characters), shifting all subsequent columns left by 8
characters.  Angle precision also differed (`{:6.2f}` vs `{:6.3f}`).  The
`[Dimer]` header was also missing a format specifier for the "Align" column.
