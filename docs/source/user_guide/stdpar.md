# nvc++ stdpar (NEB image forces)

`OpenCC` in the SURF notes is **nvc++**, not a library. GPU work there is
`nvc++ -mp=gpu` / `-stdpar=gpu`, and for QMCPACK `QMC_GPU="openmp;cuda"`.
That last flag is a QMCPACK CMake option. It is not an eOn meson option.

## What eOn offloads

The NEB image-force loop is already `std::for_each(std::execution::par, …)`
behind `-DEON_PARALLEL_NEB`. Two ways to turn that on:

| Meson option | Compiler | Backend |
|---|---|---|
| `-Dwith_parallel_neb=true` | GCC / Clang | TBB (`libstdc++`) |
| `-Dstdpar=cpu` | nvc++ | `-stdpar=multicore` |
| `-Dstdpar=gpu` | nvc++ | `-stdpar=gpu` plus `-gpu=` |

`-Dstdpar=gpu` requires Meson compiler id `nvidia_hpc` or `pgi`. GCC/Clang
configure fails rather than passing a flag the linker will reject.

```bash
CXX=nvc++ meson setup /tmp/eon-stdpar \
  -Dstdpar=gpu -Dstdpar_gpu_cc=cc80,cc90 \
  -Db_pie=false -Ddefault_library=static
```

nvc++ 23.7 does not implement Meson's `b_pie`. Leave it off. Meson's
PGICompiler has `get_pic_args` and no `get_pie_args`, so Highway's
`hwy_list_targets` dies writing ninja. eOn skips the cmake Highway wrap
on `nvidia_hpc`/`pgi`. The PGI linker also has no `link_whole`, so the
configure must be `-Ddefault_library=static`. The login-node probe needs
`libatomic.so.1` from GCCcore on `LD_LIBRARY_PATH`. Put the build dir on
a local disk: NFS home on Elja stamps files ~100 s in the future and
Meson then refuses `coredata.dat`.

The default GPU arch list is `cc80,cc90` (A100 and H100, the Snellius/Elja
pair). Override with `-Dstdpar_gpu_cc=cc90` or `-Dstdpar_gpu_cc=native` to
omit `-gpu=`. That is the same lesson as QMCPACK `QMC_GPU_ARCHS=sm_80;sm_90`,
written in nvc++'s flag language.

`parallel = true` in `config.ini` still fans out with one `std::thread` per
image when `EON_PARALLEL_NEB` is off.

## What stays on the host

The potential is the cost. Morse, LJ, and the other in-tree host potentials
keep their atom loops on the CPU. Wrapping those loops in
`std::for_each(std::execution::par)` is not the surf-notes path: nvc++
stdpar wants a larger data-parallel kernel (a whole image force, or a
GPU potential), not a three-coordinate inner loop.

`gprd_linalg_backend=stdpar` is a separate GPR-dimer linear-algebra
choice. It does not turn on NEB image offload.

## QMCPACK vs eOn

On Snellius/Elja the QMCPACK GPU build is `QMC_GPU="openmp;cuda"` with
nvc++ and CUDA 12. That string does not belong in eOn's meson options.
eOn has no OpenMP target offload and no `QMC_GPU` equivalent.

The SURF notes also record that device bitcode on disk is not evidence of
offload (`llvm-config --targets-built` must list NVPTX) and that compute
nodes without `glibc-devel` need `--sysroot` for a host compile. Those are
site toolchain facts. They do not change eOn's meson options.

Build the nvc++ tree on a login node that has headers. Do not compile on a
compute node that lacks `features.h`.

On Elja the login-node probe also needs `libatomic.so.1` from GCCcore on
`LD_LIBRARY_PATH`. `nvc++` 23.7 identifies to Meson as `nvidia_hpc`. Without
that library, Meson reports that nvc++ executables are not runnable.
