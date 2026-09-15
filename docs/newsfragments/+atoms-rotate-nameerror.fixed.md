`rotm`, `get_rotation_matrix`, and `internal_motion` call `numpy.cos` /
`sin` / `arccos`. They used to call bare `cos`/`sin`/`acos` and raise
NameError. `get_mappings` uses vector PBC distances.
