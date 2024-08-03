# Contributing

We welcome all kinds of contributions, from examples to documentation and code.
Please ensure that the testing, profiling, or leak handling instructions are
followed to facilitate a fast and efficient review process.

## Testing

Our CI runs the test suite in `debug` mode. This is important, since at runtime
many custom error messages are dropped.

```{code-block} bash
meson setup bbdir --prefix=$CONDA_PREFIX --libdir=lib -Dbuildtype=debug
meson test -C bbdir
```

## Profiling

For new contributions or refactors, test against either existing `main` (target branch) or / and `svn` `eOn`:

```{code-block} bash
# For a pretty screenshot
hyperfine -N "eonclient -s client/example/pos.con -p lj" "~/SVN/eonSVNonly/client/eonclient -s client/example/pos.con -p lj"
```

Additionally, the PR should include the tabulated result:

```{code-block} bash
hyperfine -N "eonclient -s client/example/pos.con -p lj" "~/SVN/eonSVNonly/client/eonclient -s client/example/pos.con -p lj" --export-markdown tmp
```

e.g. this [structure comparison refactor](https://github.com/TheochemUI/eOn/pull/150).

## Leak chesk

At the very least, `valgrind` (with `helgrind`) should be run on prospective
changes to ensure no leaks are introduced.

```{code-block} bash
valgrind --log-file="val_lj_sp" --leak-check=full --show-leak-kinds=all --track-origins=yes  eonclient -s client/example/pos.con -p lj
valgrind --log-file="hel_lj_sp" --tool='helgrind' eonclient -s client/example/pos.con -p lj
cat val_lj_sp
cat hel_lj_sp
```

The results should be clean, i.e. no errors (suppressed or otherwise).

```{code-block} bash
==564195== Memcheck, a memory error detector
==564195== Copyright (C) 2002-2024, and GNU GPL'd, by Julian Seward et al.
==564195== Using Valgrind-3.23.0 and LibVEX; rerun with -h for copyright info
==564195== Command: eonclient -s client/example/pos.con -p lj
==564195== Parent PID: 4149992
==564195==
==564195==
==564195== HEAP SUMMARY:
==564195==     in use at exit: 0 bytes in 0 blocks
==564195==   total heap usage: 4,085 allocs, 4,085 frees, 272,727 bytes allocated
==564195==
==564195== All heap blocks were freed -- no leaks are possible
==564195==
==564195== For lists of detected and suppressed errors, rerun with: -s
==564195== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)
```

Note that for helgrind, since we use `spdlog` one false positive [is expected](https://github.com/gabime/spdlog/issues/1425).

```{code-block} bash
==577676== Helgrind, a thread error detector
==577676== Copyright (C) 2007-2024, and GNU GPL'd, by OpenWorks LLP et al.
==577676== Using Valgrind-3.23.0 and LibVEX; rerun with -h for copyright info
==577676== Command: eonclient -s client/example/pos.con -p lj
==577676== Parent PID: 4149992
==577676==
==577676== ---Thread-Announcement------------------------------------------
==577676==
==577676== Thread #1 is the program's root thread
==577676==
==577676== ----------------------------------------------------------------
==577676==
==577676== Thread #1: pthread_cond_{signal,broadcast}: dubious: associated lock is not held by any thread
==577676==    at 0x484F1BD: pthread_cond_signal_WRK (hg_intercepts.c:1567)
==577676==    by 0x4B2C0AA: spdlog::details::periodic_worker::~periodic_worker() (in /home/rgoswami/micromamba/envs/eongit/lib/libspdlog.so.1.13.0)
==577676==    by 0x4B097A1: spdlog::details::registry::~registry() (in /home/rgoswami/micromamba/envs/eongit/lib/libspdlog.so.1.13.0)
==577676==    by 0x4CC2FA0: __run_exit_handlers (exit.c:108)
==577676==    by 0x4CC306D: exit (exit.c:138)
==577676==    by 0x4CA9C8E: (below main) (libc_start_call_main.h:74)
==577676==
==577676==
==577676== Use --history-level=approx or =none to gain increased speed, at
==577676== the cost of reduced accuracy of conflicting-access information
==577676== For lists of detected and suppressed errors, rerun with: -s
==577676== ERROR SUMMARY: 1 errors from 1 contexts (suppressed: 11 from 9)
```
