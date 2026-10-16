# H2 point job (OpenCPMD)

An ASE `molecule("H2")` in a 7.26 Å box (`pos.con`, vacuum 3.63 Å) and the
public Troullier–Martins BLYP pseudopotential `pseudo/H_MT_BLYP.psp`.

The pseudopotential is `pseudo/H_MT_BLYP.psp` from
[OpenCPMD/pseudopotentials-extended-library](https://github.com/OpenCPMD/pseudopotentials-extended-library)
(`pseudo.tar.gz`). It is norm-conserving, built for BLYP, and the library
lists `LMAX=S LOC=S` with a minimum cutoff of 25 Ry.

`scripts/ci/opencpmd_point.sh` builds pinned OpenCPMD and `libcpmdc.so`, then
runs this geometry as an eOn `point` job through
`examples/akmc-cpmd-slurm/potfiles/cpmd_extpot.py`.
