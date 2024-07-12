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

## Second approximation

The eventual goal is to have complete locality of options with their classes.
This means we want to define the parameters and defaults as close to the objects
as possible. To ensure this we have each concrete potential define a `Params`
public struct and use that instead.

### Example: LJ Parameters

Consider the iterative improvement for the `LJ` class. First the parameters need to become part of the `Parameters` object.

```{code-block} cpp
// [Potential] //
struct Potential {
  PotType potential;
  ...
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

Naturally post C++11/C++14 we can use in-class member initialization, aka
[NSDMI](https://www.cppstories.com/2015/02/non-static-data-members-initialization/).

```{code-block} cpp
namespace eonc::def {
struct LJParams {
  double u0 = 1.0;
  double cutoff = 15.0;
  double psi = 1.0;
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

The final improvement discussed earlier is to ensure complete locality, becoming:

```{code-block} cpp
class LJ : public Potential<LJ> {
public:
  struct Params {
    double u0{1.0};
    double cutoff_R{15.0};
    double psi{1.0};
    Params(const toml::table &tbl) {
      u0 = tbl["u0"].value_or(u0);
      cutoff_R = tbl["cutoff"].value_or(cutoff_R);
      psi = tbl["psi"].value_or(psi);
    }
  };

private:
  double u0;
  double cutoff_R;
  double psi;
  double cutoff_U{0};

  double calc_cutoffU(const Params &p);

public:
  LJ(const Params &ljp)
      : u0{ljp.u0},
        cutoff_R{ljp.cutoff_R},
        psi{ljp.psi},
        cutoff_U{calc_cutoffU(ljp)} {}

  void force(const ForceInput &params, ForceOut *efvdat) override final;

  void setParameters(const Params &ljp);
};
```

Where we can see that the input is always via the structure, since with
aggregate named initialization we eliminate a host of positional argument
errors.

```{code-block} cpp
case PotType::LJ: {
  // auto params = LJ::Params(config);
  // or
  auto params = LJ::Params(.u0=1.0, .cutoff_R=15.0, .psi=1.0);
  return (std::make_shared<LJ>(params));
  break;
}
```
