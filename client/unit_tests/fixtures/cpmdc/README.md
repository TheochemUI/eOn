# CPMD interface fixture

These MIT-licensed fixture sources come from rgpot commit
`e40fd175f5230098d56333bf2aff38cd839dc810`. They supply the typed CPMD
engine ABI for `test_coverage_paths`, including the existing large-box
Lennard-Jones geometry and bead-batch branches. The fixture returns
prescribed energies, forces and stresses for interface checks.

The physical CPMD engine remains selected by the environment for
`test_rgpot_pot` and the MPI backend tests. Its pseudopotentials and
scientific results are separate from this deterministic fixture.
