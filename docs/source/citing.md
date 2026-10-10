---
myst:
  html_meta:
    "description": "How to cite eOn: the Zenodo record, the method papers, and the older papers."
    "keywords": "eOn citation, Zenodo, cite eOn"
---
# Citing eOn

Please cite the eOn software through its Zenodo record. The concept DOI
[10.5281/zenodo.18529495](https://doi.org/10.5281/zenodo.18529495) always
resolves to the latest release, and each release also has its own DOI on Zenodo.
GitHub reads the same reference from `CITATION.cff`.

```bibtex
@software{goswamiEOn,
  title = {{TheochemUI/eOn}},
  author = {Rohit Goswami and Sam Chill and Rye Terrell and Graeme Henkelman and Matthew Welborn and Liang Zhang and Andreas Pedersen and {t-brink} and Erik Edelmann and Jean Claude and {Alejandro} and Sung Hoon Jung and Seyed Alireza Ghasemi and {chemist29} and {satishskamath} and {Maxim} and Tobias Dijkhuis and {via9a}},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.18529495},
  url = {https://doi.org/10.5281/zenodo.18529495}
}
```

## Papers on eOn 2.x and later

The rewrites behind eOn 2.x and 3.x, including the Gaussian process methods, are
covered in R. Goswami, *Efficient Exploration of Chemical Kinetics*, arXiv:2510.21368 (2025).
[doi:10.48550/arXiv.2510.21368](https://doi.org/10.48550/arXiv.2510.21368)

If you use one of these methods, please also cite the paper behind it.

- GPR accelerated dimer (`min_mode_method = gprdimer`, `gp_surrogate`):
  R. Goswami, M. Masterov, S. Kamath, A. Peña-Torres and H. Jónsson,
  *J. Chem. Theory Comput.* **21**, 7935 (2025).
  [doi:10.1021/acs.jctc.5c00866](https://doi.org/10.1021/acs.jctc.5c00866)
- Pruning and optimal transport for GP saddle searches (`use_prune` and
  related keys): R. Goswami and H. Jónsson, *ChemPhysChem* **27**, e202500730
  (2026). [doi:10.1002/cphc.202500730](https://doi.org/10.1002/cphc.202500730)
- Off-path climbing image NEB (OCI-NEB) with Hessian eigenmode alignment:
  R. Goswami, M. Gunde and H. Jónsson, *Front. Chem.* **14**, 1807063 (2026).
  [doi:10.3389/fchem.2026.1807063](https://doi.org/10.3389/fchem.2026.1807063)
- Benchmarking saddle searches, and rotation removal in the dimer: R. Goswami,
  *AIP Adv.* **15**, 085210 (2025). [doi:10.1063/5.0283639](https://doi.org/10.1063/5.0283639)
- Background on GP acceleration of saddle and minimum searches: R. Goswami,
  *ACS Phys. Chem. Au* **6**, 633 (2026). [doi:10.1021/acsphyschemau.6c00038](https://doi.org/10.1021/acsphyschemau.6c00038)
- NEB workflows with eOn and machine-learned potentials: R. Goswami,
  *MethodsX* **16**, 103899 (2026). [doi:10.1016/j.mex.2026.103899](https://doi.org/10.1016/j.mex.2026.103899)
- Two-dimensional RMSD plots of NEB paths (via `rgpycrumbs`): R. Goswami,
  *MethodsX* **16**, 103851 (2026). [doi:10.1016/j.mex.2026.103851](https://doi.org/10.1016/j.mex.2026.103851)

## Older papers

These describe the original EON code and the methods it started from, oldest first.

- G. Henkelman and H. Jónsson, "Long time scale kinetic Monte Carlo
  simulations without lattice approximation and predefined event table",
  *J. Chem. Phys.* **115**, 9657 (2001). [doi:10.1063/1.1415500](https://doi.org/10.1063/1.1415500)
- L. Xu and G. Henkelman, "Adaptive kinetic Monte Carlo for first-principles
  accelerated dynamics", *J. Chem. Phys.* **129**, 114104 (2008).
  [doi:10.1063/1.2976010](https://doi.org/10.1063/1.2976010)
- A. Pedersen, G. Henkelman, J. Schiøtz and H. Jónsson, "Long time scale
  simulation of a grain boundary in copper", *New J. Phys.* **11**, 073034 (2009).
  [doi:10.1088/1367-2630/11/7/073034](https://doi.org/10.1088/1367-2630/11/7/073034)
- A. Pedersen and H. Jónsson, "Distributed implementation of the adaptive
  kinetic Monte Carlo method", *Math. Comput. Simul.* **80**, 1487 (2010).
  [doi:10.1016/j.matcom.2009.02.010](https://doi.org/10.1016/j.matcom.2009.02.010)
- S. T. Chill, M. Welborn, R. Terrell, L. Zhang, J.-C. Berthet, A. Pedersen,
  H. Jónsson and G. Henkelman, "EON: software for long time simulations of
  atomic scale systems", *Modelling Simul. Mater. Sci. Eng.* **22**, 055002
  (2014). [doi:10.1088/0965-0393/22/5/055002](https://doi.org/10.1088/0965-0393/22/5/055002)
