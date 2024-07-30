# Tests

```{versionadded} 2.0
```

We have a robust testsuite, consisting of unit tests, approval tests, and a few
integration tests.

## Writing and Registering Tests

We find that, rather than build each executable by hand or even register each
one by hand, we can leverage the array iteration features of the `meson`
language[^1].

```{code-block} meson
if get_option('build_tests')
  _args += ['-DEONTEST'] # Unused
  _deps += [ gtest_dep ]
test_array = [#
  ['Improved Dimer', 'impldim_run', 'ImpDimerTest.cpp', '/gtests/data/saddle_search'],
             ]
foreach test : test_array
  test(test.get(0),
       executable(test.get(1),
          sources : ['gtests/'+test.get(2)],
          dependencies : [ eon_deps, gtest_dep, gmock_dep ],
          link_with : eclib,
          cpp_args : eon_extra_args
                 ),
        workdir : meson.source_root() + test.get(3)
      )
endforeach
endif
```

[^1]: Described [here](https://click.rgoswami.me/auto-disc-meson-tests)
