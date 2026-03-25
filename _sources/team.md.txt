---
myst:
  html_meta:
    "description": "Details of the eOn software development team and key contributors."
    "keywords": "eOn team, developers, contributors"
---
# eOn Team

`eOn` is maintained by [Rohit Goswami](https://rgoswami.me) and the
[Jónsson Group](https://hj.hi.is/researchgroup.html) at the University of
Iceland.  For background on the project's history and relationship to the
UT Austin codebase, see
[this post](https://rgoswami.me/posts/eon-acad-foss/).

## Active Contributors

### [Jónsson Group](https://hj.hi.is/researchgroup.html) (University of Iceland)

- [Rohit Goswami](https://rgoswami.me), **maintainer & release manager**
- Andreas Pederson
- Alejandro Pena Torres
- Vilhjálmur Ásgeirsson

## Historical Contributors (pre-v2.8)

The Henkelman group were major contributors up to v2.8.  See
[their project page](https://theory.cm.utexas.edu/eon/) and
[this post](https://rgoswami.me/posts/eon-acad-foss/) for background.

### [Henkelman Group](https://theory.cm.utexas.edu/henkelman/) (UT Austin)

- Samuel Chill
- Rye Terrell
- Matthew Welborn
- Liang Zhang
- Akksay Singh
- S. Alireza Ghasemi
- Juliana Duncan
- Naman Katyal
- Sung Hoon Jung

### [Li Lei Group](https://faculty.sustech.edu.cn/lil33/en/)

- SONG Zichen

### Other contributions

With contributions from a number of other developers:

- Tobias Brink (TU Darmstadt)
- Ian Johnson (UT)
- Saif Kazim (UT)
- Jean Claude Berthet (Iceland)
- Raymond Smith (UT)
- Erik Edelmann (Finland)
- Jutta Rogal (Bochum)

### Related software and collaborators

The development of eOn has been influenced and supported by the following key
collaborators and groups, many of which came together for the
[CECAM-LTS-MAP](https://github.com/CECAM-LTS-MAP) workshop:

- The [MAMASMIAS consortium](https://mammasmias.gitlab.io/), represented by [Dr.
  Miha Gunde](https://github.com/mgoonde), has had a major impact on recent
  versions of the code, particularly in improving its interoperability.[^1]

Discussions on the software's architecture have been advanced by our
collaborators in Canada:

- [Prof. Laurent Béland](https://www.laurentkarimbeland.com/) (Queen's
  University, Kingston) led strategic discussions on the design of KMC
  algorithms and provided direct support for the project by funding Rohit
  Goswami's development work for a quarter in 2025 towards combining
  disparate KMC codes[^2] .

    - The k-ART[^3]  group of [Prof. Normand
      Mousseau](https://phys.umontreal.ca/repertoire-departement/professeurs/professeur/in/in14441/sg/Normand%20Mousseau/)
      (University of Montreal) contributed to discussions on project's
      high-level design.

Rohit also acknowledges funding support from
[lab-COSMO](https://www.epfl.ch/labs/cosmo/), for supporting
[metatomic](https://docs.metatensor.org/metatomic/latest/overview.html) models
in eOn, in particular [Dr. Guillaume Fraux](https://guillaume.fraux.fr/) and
[Prof. Michele Ceriotti](https://people.epfl.ch/michele.ceriotti?lang=en).

## License history

`eOn` was originally licensed under the GNU GPL v3 until 2017.  In 2017, the
license was changed to GNU GPL v2.  Since 2020, eOn is licensed under the BSD
3-Clause License. Some more context is at
[gh-140](https://github.com/TheochemUI/eOn/issues/140).

### Footnotes

[^1]: [IRA](https://mammasmias.github.io/IterativeRotationsAssignments/intro.html) has also been used in tandem with eOn for research in 2024-2025

[^2]: The work-in-progress (2025) pyKMC code blends concepts from eOn and k-ART

[^3]: [k-ART](https://kart-doc.readthedocs.io/en/latest/download.html) is an off-lattice KMC code primarily used for metals but with overlapping concepts
