``displace_atom_list`` is CON file-order. ``ListedAtoms`` always remaps
those rows through the ``atom_id`` sort, then keeps free atoms, so a
movable-first active-volume ``.con`` no longer raises "Listed atoms
are all frozen" or displaces a coincidentally-free Structure row.
