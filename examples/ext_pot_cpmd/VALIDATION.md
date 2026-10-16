# External CPMD example validation

The point and NEB commands were exercised with a real OpenCPMD executable,
not a substituted potential. The client was conda-forge eOn 3.6.0
(`c4efc62`); OpenCPMD was `062582b7cfd832d36f88f504cd08e4ead42eb404`,
built with GNU Fortran 14.4 and Open MPI 5.0.11. The example scripts were
checked at upstream commit `95f934464`.

## Single point from CON and XYZ

The supplied seven-silicon vacancy structure, evaluated with neutral BLYP,
70 Ry and `Si_MT_BLYP.psp` with `LMAX=D LOC=D`, gives
**−739.9932635496 eV**. Both the original CON and an extended XYZ carrying
the same cell give this value through `run.sh`. The XYZ path invokes the
pinned `readcon-chemfiles==0.17.1` converter through uv.

This reference belongs to that exact input and pseudopotential. A different
structure, cutoff or potential file has a different reference energy.

## NEB through the same executable

A periodic translation provides a check whose expected barrier is zero:
the product translates all seven atoms by 0.125 Angstrom along each axis,
without changing their relative positions or cell. The test uses BLYP,
20 Ry and `Si_MT_BLYP.psp` with `LMAX=P LOC=P`, matching the supplied
vacancy's relaxation settings.

`run_neb.sh` with the distributed five-image configuration reports:

```text
GOOD termination_reason_text
7 total_force_calls
5 number_of_images
```

The largest internal-image force is **0.008902 eV/Angstrom**, below the
configured 0.05 eV/Angstrom criterion. Image energies span approximately
**0.000076 eV**, consistent with a numerical check of a zero barrier.
Both CON endpoints and extended-XYZ endpoints were exercised. This checks
the working NEB/oracle interface; it does not measure a reactive saddle.

## File protocol and conversion checks

A separate prescribed-output fixture checks the wrapper independently:
−2 Hartree becomes −54.42277249198 eV, forces are multiplied by
51.4220674763, and interleaved Si/N atoms return in their original order.
This fixture is a protocol test, not a DFT calculation. eOn can remove the
common net force, so its printed force norm is compared with that projection;
the wrapper's raw force output is also checked separately.

Readcon round trips cover multi-frame selection, Angstrom/Bohr input,
cell metadata and explicit overrides, loud warnings for inferred cells,
and preserved fixed-atom flags. NEB preparation checks a shared cell and
translation and rejects mismatched cells or atom identities.

The Intel oneAPI module names and launcher in `.env.example` are site
settings. Their execution must be checked on that site. `cpmd.out` reports
launcher initialization, module setup and executable launch, so a site
failure can be distinguished from an SCF or input failure.
