---
myst:
  html_meta:
    "description": "How to add an empirical potential to eOn by landing the kernel in rgpot."
    "keywords": "eOn, rgpot, potential, RgpotAdapter, porting"
---

# Porting potentials

Empirical kernels live in [rgpot](https://github.com/OmniPotentRPC/rgpot).
eOn does not grow an in-tree `NewPot`. A kernel that belongs next to LJ,
Morse, or the Fortran pots is added in rgpot, then selected from eOn.

## Kernel

1. Add the potential under rgpot (`CppCore/rgpot/<Name>/`) and link it
   into the single `librgpot` SONAME. Fortran kernels use the Fortran
   2018 layout already in that tree and keep their symbols hidden.
2. Expose a C++ face (`force` / `energy`) and a small config struct.
   Neighbour finding goes through rgpot's vesin cache, not a private
   pair loop.
3. Pin a reference energy and force in an rgpot test.

## eOn arm

1. Map a `potential` token to that type in `makePotential`
   (`client/Potential.cpp`) with `makeRgpot` or `makeRgpotDefault`.
2. Add the token to the parameter schema (`[Potential]` / the pot's own
   section) and to `PotType`.
3. Keep subprocess and file engines (VASP, LAMMPS, ASE, AMS, socket
   NWChem) in eOn. They are not kernels.

Direct NWChem, CPMD, metatomic, and xTB engines stay `dlopen` plugins
behind potential type `RGPOT`. Non-Windows builds always link that arm.
There is no `-Dwith_rgpot` switch.

`EON_POTENTIALS_PATH` and `[Potential] potentials_path` still name
directories searched for those engine plugins. They are not a way to
register a new empirical kernel.
