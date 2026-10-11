A convel file now labels velocities with the natural time unit (amu*angstrom^2/eV)^0.5 and writes con_spec_minor 1.
eOn 3.6.0 and earlier (doi:10.5281/zenodo.23166734) stored those velocities under units.time fs and wrote no minor.
eOn still reads spec versions 1, 2, and 3.0 as the natural unit, so those files stay correct.
readcon 0.17 does not apply that fs label. readcon-db export_h5md takes velocity_time_unit set to (amu*angstrom^2/eV)^0.5 when a later file mislabels the same numbers as fs.
Convel text is binary64 round-trip unless a decimal count is requested. Six decimal places lose about 5.8e-7 relative kinetic energy.
