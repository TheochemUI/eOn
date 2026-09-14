``displace_atom_list`` is CON file-order. ``ListedAtoms`` and the C++
``listedAtomEpiCenter`` path always remap those rows through the
``atom_id`` sort, then keep free atoms. Lone ``-1`` is every free atom.
A movable-first active-volume ``.con`` no longer raises "Listed atoms
are all frozen" or displaces a coincidentally-free post-sort row.
Indices from ``displace_atom_kmc_state_script`` stay Structure rows of
the temp ``savecon`` file and are not remapped as original file-order.
