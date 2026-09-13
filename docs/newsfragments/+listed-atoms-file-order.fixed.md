``displace_atom_list`` is CON file-order. ``ListedAtoms`` and the C++
``listedAtomEpiCenter`` path always remap those rows through the
``atom_id`` sort, then keep free atoms. Lone ``-1`` is every free atom.
A movable-first active-volume ``.con`` no longer raises "Listed atoms
are all frozen" or displaces a coincidentally-free post-sort row.
