``listedAtomEpiCenter`` now picks with ``randomDouble(size)`` so the last
listed (or last free, for lone ``-1``) atom can be chosen. The previous
``size-1`` interval always dropped that last index.
