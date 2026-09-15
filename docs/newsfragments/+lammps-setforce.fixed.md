LAMMPS applies `fix setforce 0 0 0` to atoms with `Atom.fixed` set on all three axes before `run 1`, so buffer atoms stay put during post-saddle minimization.
