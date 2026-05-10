---
myst:
  html_meta:
    "description": "Guide for developers on the testing suite in eOn, including unit tests, approval tests, and integration tests, and how to add new tests."
    "keywords": "eOn testing, unit tests, integration tests, meson tests, Catch2"
---

# Tests

```{versionchanged} 2.12
Migrated from GoogleTest to Catch2. All tests use the `TEST_CASE` /
`TEST_CASE_METHOD` macros with Catch2 matchers.
```

We have a robust test suite consisting of unit tests and integration tests.
All tests are registered via meson's array iteration pattern and run with
`meson test -C builddir`.

## Test inventory

| Test | File | Data | What it covers |
|---|---|---|---|
| `test_impldim` | `ImpDimerTest.cpp` | `saddle_search` | Improved dimer eigenmode search |
| `strparse_run` | `StringHelpersTest.cpp` | -- | String parsing utilities |
| `test_matter` | `MatterTest.cpp` | `neb_morse` | Matter construction, copy, PBC, forces |
| `test_pot` | `PotTest.cpp` | `neb_morse` | LJ potential energy/force sanity |
| `test_zbl` | `ZBLPotTest.cpp` | `au_si` | ZBL potential switching function |
| `test_extpot` | `ExtPotTest.cpp` | `extpot` | External potential file interface |
| `test_proj_rot_trans` | `ProjectOutRotTransTest.cpp` | -- | Rotation/translation projection |
| `test_epicenters` | `EpiCentersTest.cpp` | `Pt_Heptamer_FrozenLayers` | Displacement epicenter selection |
| `test_storage_order` | `StorageOrderTest.cpp` | -- | Eigen storage order verification |
| `test_pot_registry` | `PotRegistryTest.cpp` | -- | Force call tracking registry |
| `test_neb_lj` | `NEBMorseTest.cpp` | `neb_morse` | NEB convergence on Morse cluster |
| `test_neb_regression` | `NEBRegressionTest.cpp` | `neb_morse` | NEB regression (CI-NEB barrier) |
| `test_dimer` | `DimerTest.cpp` | `neb_morse` | Eigenmode strategy pattern |
| `test_saddle` | `SaddleSearchTest.cpp` | `neb_morse` | MinModeSaddleSearch integration |
| `test_confileio` | `ConFileIOTest.cpp` | `neb_morse` | CON file round-trip I/O |
| `test_socket_nwchem` | `SocketNWChemPotTest.cpp` | `nwchem_test` | NWChem socket potential (needs nwchem) |

Optional tests (enabled by build flags):
- `test_ase_pot` (with_ase): ASE Python calculator
- `test_mta` (with_metatomic): Metatomic ML potential
- `test_serve_spec` (with_serve): Serve mode endpoint parsing

Runtime-loaded tests (always built; SKIP at runtime when the .so is
absent from `LD_LIBRARY_PATH`):
- `test_xtb`, `test_cineb_xtb`: libxtb
- `test_artn`: libartn
- `test_ira`: libira

## Writing and registering tests

Tests use [Catch2](https://github.com/catchorg/Catch2) (amalgamated, vendored
in `thirdparty/catch2/`). Each test file is compiled with the amalgamated
source and linked against `eonclib`.

Registration uses meson's array iteration:

```{code-block} bash
:caption: client/meson.build (excerpt)
test_array = [
    ['test_name', 'TestFile.cpp', 'data_dir'],
]
foreach test : test_array
    test(test.get(0),
        executable(test.get(0),
            sources: ['gtests/' + test.get(1),
                      'thirdparty/catch2/catch_amalgamated.cpp'],
            dependencies: test_deps,
            include_directories: _incdirs,
            cpp_args: test_args,
            link_with: _linkto,
        ),
        workdir: meson.project_source_root()
            + '/client/unit_tests/data/systems/' + test.get(2),
    )
endforeach
```

## Test data

Test data lives in `client/unit_tests/data/systems/`. Each subdirectory contains
the `config.ini`, `pos.con`, and any other files needed by the test. The
`neb_morse` dataset (7-atom LJ cluster) is shared by most tests.

## Running tests

```bash
# All tests
pixi run -e dev meson test -C builddir

# Single test
pixi run -e dev meson test -C builddir test_neb_lj

# With verbose output
pixi run -e dev meson test -C builddir -v
```
