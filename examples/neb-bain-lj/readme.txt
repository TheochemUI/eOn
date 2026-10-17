Bain path of a Lennard-Jones crystal: bcc to fcc by one tetragonal strain,
on the solid-state band (solid_state = true). Run with eonclient in this
directory.

The potential is LJ with u0 = 1 eV and psi = 1 Angstrom (the species label
Ar only names the atoms; these are not argon's parameters), cut at 2.2
Angstrom with the C^2 switch over the last 0.2 Angstrom. Under that
truncation bcc is a local minimum; with a cutoff of 2.5 Angstrom or more
it is not. The endpoints are 5 x 5 x 4 two-atom body-centred tetragonal
cells (200 atoms), written by make_endpoints.py. Every atom sits on an
inversion centre, so the path is the cell strain alone.

Reference, from a relaxed scan of the energy per atom over c/a with spacing
2.5e-4 (validation/bain_lj_reference.py in rgsaddle):

  bcc     c/a 1           E = -7.4148256256 eV/atom
  saddle  c/a 1.09331985  E = -7.3557899883 eV/atom
  fcc     c/a sqrt(2)     E = -7.6698410342 eV/atom
  barrier 0.0590356373 eV/atom, 11.80712746 eV for 200 atoms

The climbing image should reach c/a = (h22 / 4) / (h00 / 5) within 2.5e-4
of 1.09332, and a barrier within 1.2e-6 eV/atom of the scan, the energy
change over one scan step at the saddle. client/unit_tests/
SolidStateBainTest.cpp runs the same band and checks both.
