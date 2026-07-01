`ConFileIO` targets readcon-core **0.13**: per-axis fixed masks, optional forces
and frame energy when the potential cache is current, `metadata_from_frame`
(NEB / potential_type / energy fields), `z_to_symbol`, and compressed writers via
`ConFrameWriter::compression_from_extension` (`.con.gz` / `.con.zst`). Subproject
wrap pins `v0.13.1`.
