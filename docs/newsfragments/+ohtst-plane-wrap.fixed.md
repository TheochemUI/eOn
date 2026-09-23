OH-TST keeps the unwrapped free coordinates in `samplePlane` when the next hyperplane is loaded. A reload after `setPositionsFreeV` shifted a coordinate by a lattice vector, and the projection then left the constrained restart.

