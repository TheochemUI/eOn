# QUIP with LAMMPS

```{versionadded} 2.0
```

In order, we require:
- QUIP as a library (`quippy` is not required)
  + SOAP needs to be compiled in, [documented here](https://github.com/HaoZeke/quip-nix)
- LAMMPS as a shared library, [documented here](https://docs.lammps.org/Build_link.html)
  + Without MPI
  + Linked to the same QUIP
  
The rest of the document walks through a concrete example of this use-case.

## Sample run

This needs to be in conjunction with the LAMMPS [pair style
quip](https://docs.lammps.org/pair_quip.html).

## Example Build
  
## Base environment

```{code-block} bash
micromamba create -n quip_eon_lammps
micromamba activate quip_eon_lammps
micromamba install -c conda-forge eigen numpy ase PyYAML pytest sh pytest-datadir sphinx meson gh cmake subversion pkgconfig openblas
export PYTHONROOT=$CONDA_PREFIX
```

## QUIP

```{code-block} bash
git clone https://github.com/HaoZeke/quip-nix.git
cd quip-nix
git clone --recursive https://github.com/libAtoms/QUIP.git QUIP
cd QUIP/src
git clone --recursive https://github.com/mcaroba/turbogap.git 
source ./../setEnvVars.sh
mkdir -p build/$QUIP_ARCH 
mkdir -p "$QUIP_STRUCTS_DIR"
cp ../files/$QUIP_ARCH-Makefile.inc build/$QUIP_ARCH/Makefile.inc 
export EXTRA_LINKOPTS="$(pkg-config --libs openblas)"
make -j1
make libquip
make install 
```

We need to also keep track of the build system outputs as well (though it is
also installed into the environment in theory by now).

```{code-block} bash
cd build/linux_x86_64_gfortran/
export QUIP_LIB_DIR_BLD=$(pwd)
```

## LAMMPS

We require serial `LAMMPS`, linked to our recently installed `QUIP`.

```{code-block} bash
git clone https://github.com/lammps/lammps
cd lammps
mkdir build; cd build
# The QUIP_LIB_DIR_BLD is special, needs to be populated from the QUIP folder
cmake ../cmake -D CMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DDOWNLOAD_QUIP=no -DQUIP_LIBRARY=$QUIP_LIB_DIR_BLD/libquip.a -DBUILD_SHARED_LIBS=yes -D PKG_ML-QUIP=yes -D BUILD_MPI=no -D PKG_MANYBODY=yes
make -j$(nproc)
make install
```

Assuming the same environment is not being used, it makes sense to track the
variables for the outputs as well.

```{code-block} bash
cd build
export LMP_BLD_DIR=$(pwd)
```

## eOn

Assuming we are starting from scratch..

```{code-block} bash
gh repo clone theochemui/eon 
cd eOn
ln -sf $LMP_BLD_DIR/liblammps* "./client/potentials/LAMMPS/"
meson setup bbdir --prefix=$CONDA_PREFIX --libdir=lib -Dwith_lammps=True --buildtype=release
meson install -C bbdir
```

At this point both `eonclient` examples and `eon` (e.g. AKMC) examples will run.

