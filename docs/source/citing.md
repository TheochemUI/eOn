---
myst:
  html_meta:
    "description": "How to cite eOn: the Zenodo record, the method papers, and the older papers."
    "keywords": "eOn citation, Zenodo, cite eOn"
---
# Citing eOn

Please cite the eOn software through its Zenodo record
{cite:p}`cite-goswamiTheochemUIEOn2026`. The concept DOI
[10.5281/zenodo.18529495](https://doi.org/10.5281/zenodo.18529495) always
resolves to the latest release, and each release also has its own DOI on Zenodo.
GitHub reads the same reference from `CITATION.cff`. The entry below matches
`docs/source/bibtex/eonDocs.bib`, which these pages draw on.

```bibtex
@misc{goswamiTheochemUIEOn2026,
  title = {{{TheochemUI/eOn}}},
  author = {Goswami, Rohit and Chill, Sam and Terrell, Rye and Henkelman, Graeme and Welborn, Matthew and Zhang, Liang and Pedersen, Andreas and {t-brink} and Edelmann, Erik and Claude, Jean and {Alejandro} and Jung, Sung Hoon and Ghasemi, Seyed Alireza and {chemist29} and {satishskamath} and {Maxim} and Dijkhuis, Tobias and {via9a}},
  year = 2026,
  publisher = {Zenodo},
  howpublished = {Zenodo},
  doi = {10.5281/zenodo.18529495},
  url = {https://doi.org/10.5281/zenodo.18529495},
  note = {Concept DOI, resolves to the latest release}
}
```

## Papers on eOn 2.x and later

The rewrites behind eOn 2.x and 3.x, including the Gaussian process methods and
the socket and calculator interfaces, are covered in
{cite:t}`cite-goswamiEfficientExplorationChemical2025`.

If you use one of these methods, please also cite the paper behind it.

- Gaussian process regression accelerated dimer
  (`min_mode_method = gprdimer`):
  {cite:t}`cite-goswamiEfficientImplementationGaussian2025`.
- Adaptive pruning and optimal transport for Gaussian process saddle searches
  (`use_prune` and related keys under `[GPR Dimer]`):
  {cite:t}`cite-goswamiAdaptivePruningIncreased2025b`.
- Off-path climbing image NEB (OCINEB, `ci_mmf = true`) with Hessian eigenmode
  alignment: {cite:t}`cite-goswamiEnhancedClimbingImage2026`.
- Benchmarking saddle searches, and rotation removal in the dimer
  (`remove_rotation`): {cite:t}`cite-goswamiBayesianHierarchicalModels2025a`.
- Background on Gaussian process acceleration of saddle and minimum searches:
  {cite:t}`cite-goswamiTutorialReviewBayesian2026`.
- NEB workflows with eOn and machine-learned potentials:
  {cite:t}`cite-goswamiReproducibleOrchestrationBest2026`.
- Two-dimensional RMSD plots of NEB paths:
  {cite:t}`cite-goswamiTwodimensionalRMSDProjections2026`.
- Metatomic potentials in eOn, including the NEB runs:
  {cite:t}`cite-bigiMetatensorMetatomicFoundational2026`.

## Older papers

These describe the original EON code and the methods it started from, oldest
first: kinetic Monte Carlo without a lattice
{cite:p}`cite-henkelmanLongTimeScale2001`, adaptive kinetic Monte Carlo
{cite:p}`cite-xuAdaptiveKineticMonte2008`, long time scale simulation of a
grain boundary {cite:p}`cite-pedersenLongTimeScale2009`, the distributed EON2
code {cite:p}`cite-pedersenDistributedImplementationAdaptive2010`, and the 2014
software paper {cite:p}`cite-chillEONSoftwareLong2014`.

The pages for each method, such as <project:user_guide/akmc.md>,
<project:user_guide/dimer.md> and <project:user_guide/neb.md>, cite the papers
those methods come from.

## References

```{bibliography}
---
style: unsrt
filter: docname in docnames
labelprefix: CITE_
keyprefix: cite-
---
```
