# Virtualbox Troublshooting

```{warning}
It is **highly** recommended to run `eON` locally, since hypervision virtualization incurs a **significant** slowdown.
```

The virtual machine used for the CECAM-LTS-MAP 2024 event has an older version
of `eON` which must be removed.

To do so:

```{code-block} bash
# Remove the directory
rm -rf "$HOME/software/eon"
# Remove the paths
sed -i.bak '/export PYTHONPATH=\$HOME\/software\/eon:\$PYTHONPATH/d; /export PATH=\$HOME\/software\/eon\/bin:\$HOME\/bin:\/opt\/bin:\$PATH/d' "$HOME/.bashrc"
# Refresh the shell
source "$HOME/.bashrc"
```

Finally, since most of machines are not configured with enough memory to run
`ninja` (executed by `meson`) directly, the commands should be modified to:

```{code-block} bash
meson setup bbdir --prefix=$CONDA_PREFIX --libdir=lib
# j1 means use one core only
meson compile -j1 -C bbdir
meson install bbdir
```

Now everything should be ready to run.

```{code-block} bash
cd examples/akmc-al
python -m eon.server
Eon version None
State list path does not exist; Creating: .//states/
Registering results
Processed 0 results
Queue contains 0 searches
Making 8 process searches
Job finished: .//jobs/scratch/0_0
Job finished: .//jobs/scratch/0_1
Job finished: .//jobs/scratch/0_2
```
