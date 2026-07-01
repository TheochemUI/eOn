`ConFileIO` targets readcon-core **0.13** with a nanobind-ready `IoStatus` enum
(replacing `bool` returns), bulk/zero-copy builder maps
(`positions_data` / `set_*_from_flat` / `set_atom_velocity`), and
`writeNebPath` using `ConFrameBuilder::clone()` so NEB bands write in one pass
without re-reading the movie file per image. Matter exposes atom-index and CON
header accessors for bindings.
