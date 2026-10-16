# CPMD interface fixture

These MIT-licensed fixture sources come from rgpot commit
`e40fd175f5230098d56333bf2aff38cd839dc810`. They supply the typed CPMD
engine ABI for `test_coverage_paths` and `test_rgpot_pot`. The fixture
returns prescribed energies, forces and stresses for interface checks.
`test_rgpot_pot` checks the prescribed energy `0.75 + 0.001 * cell_zz`.

The large-box Lennard-Jones geometry and bead-batch coverage branches
execute against this deterministic interface. The MPI backend tests
retain the engine selected by the calling environment. Physical CPMD
energies and forces require independent validation with their actual
method, structure and pseudopotentials.
