`ConFileIO` targets readcon-core **0.13** with the mutation-oriented builder API:
seed atoms then `set_positions_from_flat` / `set_forces_from_flat` /
`set_atom_velocity`, plus `metadata_from_frame`, `fixed_mask()`, `z_to_symbol`,
and compression-aware writers. Subproject wrap pins `v0.13.1`.
