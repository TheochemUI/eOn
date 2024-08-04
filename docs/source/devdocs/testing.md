# Tests

```{versionadded} 2.0
```

We have a robust testsuite, consisting of unit tests, approval tests, and a few
integration tests.

## Writing Approval Tests

As different functions (reading and writing) get spun out of `Matter`, the following testing protocol is to be enforced.

A simple file reader is to drive the entire equality process.

```{code-block} cpp
std::string readFileContent(const std::string &filename) {
  std::ifstream file(filename);
  std::stringstream buffer;
  buffer << file.rdbuf();
  return buffer.str();
}
``` 

Matter objects may be generated in the standard way, for writers this could be:

```{code-block} cpp
std::vector<eonc::Matter> getTestMatter() {
  const auto config =
      toml::table{{"Potential", toml::table{{"potential", "lj"}}}};
  auto pot_default = eonc::makePotential(config);
  auto matter = eonc::Matter(pot_default);
  std::string confile("pos.con");
  matter.con2matter(confile);
  return {matter};
}
```

Crucially, the **first approval** is done with the _old method_

```{code-block} cpp
TEST_CASE("VerifyMatter2Con") {
  auto testMatter = getTestMatter()[0];
  std::string filename = "test_output.con";
  testMatter.matter2con(filename);
  std::string fileContent = readFileContent(filename);
  ApprovalTests::Approvals::verify(fileContent);
}
```

Once this has been approved, then the new method / function / design is to be
used.

```{code-block} cpp
TEST_CASE("VerifyMatter2Con") {
  auto testMatter = getTestMatter();
  std::string filename = "test_output.con";
  eonc::io::ConWriter conWriter;
  conWriter.write(testMatter[0], filename);
  std::string fileContent = readFileContent(filename);
  ApprovalTests::Approvals::verify(fileContent);
}
```

In this manner, the resulting functions rigorously pass consistency checks.

## Registering Tests

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
