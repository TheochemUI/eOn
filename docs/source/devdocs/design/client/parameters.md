# Parameters

## Older design

Essentially, there was a vendored `INI` file reader, which was used to generate
a massive `Parameters` object, which was then passed around. Several comments:

- The design made it easy to not always set default values.
- The values were set multiple times, as the `INI` reader had fallbacks and the constructor also had additional initialization.
- The values were "far from their use" in that the parameters requested by each sub-part was not immediately clear.
- The object was needlessly large.


## First approximation

```{versionadded} 2.x
```

- TOML reader
- Approval tests to help with default values.
- Structures within `Parameters`

While this addressed some of the issues, it is still not clean enough.

```{todo}
An `unordered_map` is read into from the TOML file, or kept as the standard TOML
output. Each class has its own parameters, which are initialized either from the
TOML output or directly through the initializer. Since we use a factory pattern
for many objects anyway, this should be fairly straightforward.
```

## Second approximation

The eventual goal is to have complete locality of options with their classes.
This means we want to define the parameters and defaults as close to the objects
as possible.

### Example: LJ Parameters

Consider the iterative improvement for the `LJ` class. First the parameters need to become part of the `Parameters` object.

```{code-block} cpp
// [Potential] //
struct Potential {
  PotType potential;
  double MPIPollPeriod;
  bool LAMMPSLogging;
  int LAMMPSThreads;
  bool EMTRasmussen;
  bool LogPotential;
  std::string extPotPath;
  struct LJParams{
  double u0;
  double cutoff;
  double psi;
  } lj;
} pot;
```

This is already slightly better than having no scoping (as in the earlier
versions), since we don't need to have long variable names like `LJPotu0`.

However, the values are still initialized in the constructor of `Parameters`,
which is unintuitive and error-prone (since they are used in `LJ.h`).

The solution is to move the struct definition into `LJ.h`, so this now becomes:

```{code-block} cpp
// [Potential] //
struct Potential {
  PotType potential;
  ...
  eonc::def::LJParams lj;
} pot;
```

Where we also set the defaults within the `LJ.h` file itself:

```{code-block} cpp
namespace eonc::def {
struct LJParams {
  double u0;
  double cutoff;
  double psi;
  LJParams()
      : u0{1.0},
        cutoff{15.0},
        psi{1.0} {}
};
} // namespace eonc::def
```

Now we can finally remove the initialization call in the `Parameters`
constructor by addding an initializer for `Potential`.

```{code-block} cpp
// [Potential] //
struct Potential {
  PotType potential;
  ...
  eonc::def::LJParams lj;
  Potential()
      : potential{PotType::LJ},
        ...
        lj{eonc::def::LJParams()} {}
} pot;
```

The final improvement is to replace parameters completely, by a simple parser,
so as to prevent having to include headers in `Parameters.h` and quickly enable
more input file formats.

```{code-block} cpp
// TODO(rg): Remove these when parameters goes away
#include "potentials/LJ/LJ.h"
```
