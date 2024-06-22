# Installation

Eon is divided up into two separate programs: a server and a client. The client
does most of the computation (e.g. saddle searches, minimizations, and molecular
dynamics) while the server creates the input for the client and processes its
results.

# Obtaining sources

```{versionadded} 2.0
`eON` is now developed and distributed primarily via GitHub.
```

Assuming access has been granted to the Github repository:

```{code-block} bash
git clone https://github.com/TheochemUI/EONgit.git
cd EONgit
```

````{margin}
```{note}

* The [GitHub CLI tool](https://cli.github.com/) makes authentication much easier.
* Micromamba ([installation](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)) is a faster conda-helper
```
````

We provide a `conda` environment, which is only partially supported for reproducible usage, since it depends on local compilers.

```{code-block} bash
micromamba create -f environment.yml
micromamba activate eongit
```

This leads to the most robust installation approach:

```{code-block} bash
meson setup bbdir --prefix=$CONDA_PREFIX
meson install -C bbdir
```

The server is accessed through `python -m eon.server`, and the `eonclient`
binary is automatically made available in the activated environment..

```{versionchanged} 2.0
While reading older documentation, calls to `eon` must now be `python -m
eon.server`. 
```


# Additional topics

This page lists generally applicable installation instructions, for specific
systems, follow the sub-parts of this document.

```{toctree}
:maxdepth: 2
:caption: Special topics

windows
lammps
```
<!-- pipx run pdm run sphinx-build -b html docs/source docs/build/html -->

