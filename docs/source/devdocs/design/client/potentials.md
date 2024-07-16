# Potential design

Potentials in `eON` were designed for inter-operability from the start, with the business end of the implementation handed off to a generic `C-style` call:

```{code-block} cpp
AtomMatrix force(long nAtoms, AtomMatrix positions,
                 VectorXi atomicNrs, double *energy, Matrix3d box);
```

```{versionchanged} 2.x
- To support surrogate potentials, the signature include `double *variance` which is typically `nullptr` for exact potentials.
- A factory function is used for thread-safety instead of the static raw pointer `getPotential`
```

## Redesign

```{versionchanged} 3.x
```

Conceptually, the `Potential` objects are a good place for inheritance. They
have only one layer down and only one public facing function. Like the `Job`
instances, they construct parameters from `toml::table` instances directly.

### Naming

Since naming is hard, we prefer to use the following conventions:

- `Params` for the structure of parameters
  - Derived from base `PotParams` to enforce `from_toml` for setting inputs
  - Since it is a plain data object, structured designated initializers work for
    the defaults.
  - Never use `Params` without a quantifier! e.g. `LJ::Params` in arguments.
- `fip` for `ForceInput` (see below)
- `efvd` for `ForceOut` (see below)

### Parameters


```{code-block} cpp
class LJ final : public Potential<LJ> {
public:
  struct Params final : public PotParams {
    double u0{1.0};
    double cutoff_R{15.0};
    double psi{1.0};
    void from_toml(const toml::node_view<const toml::node> &tbl) override {
      u0 = tbl["u0"].value_or(u0);
      cutoff_R = tbl["cutoff"].value_or(cutoff_R);
      psi = tbl["psi"].value_or(psi);
    }
  };
};
```

Where we have to refrain from using a more elegant constructor to prevent having
to write more boilerplate (a plain data struct) for the parameters due to the
designated [initializers handling of
C++](https://stackoverflow.com/questions/64770166/why-i-can-not-use-designated-initalizers-with-structs-that-are-not-aggregates).

The `PotParams` class simply enforces the constraint of having a `from_toml` method.

```{code-block} cpp
struct PotParams {
  virtual ~PotParams() = default;
  virtual void from_toml(const toml::node_view<const toml::node> &tbl) = 0;
};
```

If there are no parameters to be set by the user, then the class doesn't have a
`Params` struct at all.

```{code-block} cpp
class Aluminum : public Potential<Aluminum> {
public:
  Aluminum() { potinit_(); };
  void forceImpl(const ForceInput &, ForceOut *) override;
};
```

### Structured force calls

The older interface was

```{code-block} c
void force(long N, const double *R, const int *atomicNrs, double *F,
           double *U, double *variance, const double *box);
```

Which has all the standard problems of having a long argument list:

- Hard to update parameters
- Are positional, but not keywords
- Make memory allocation a bit more unclear

These can elegantly be solved with structs. Additionally, designated
initializers are [part of C99](https://stackoverflow.com/a/65546688/1895378) and
[C++20](https://www.cppstories.com/2021/designated-init-cpp20/), making this a
very natural extension.

```{code-block} C
typedef struct {
  // pointer to number of atoms, pointer to array of positions
  // address to supercell size
  const size_t nAtoms;
  const double *pos;
  const size_t *atmnrs;
  const double *box;
} ForceInput;

typedef struct {
  // pointer to array of forces
  double *F;
  // Internal energy
  double energy;
  // Variance here is 0 when not needed and that's OK
  double variance;
} ForceOut;
```

One consequence of this is there is now, in this single instance, `C_Structs.h`
structures which are do not **default initialize** their values, thus setting
them apart from the structs used elsewhere, e.g. for parameters.

With this, the interface to be implemented by the class is simply:

```{code-block} C
void forceImpl(const ForceInput &, ForceOut *) override;
```

Where additional input or output parameters can be added without making
error-prone search and replace / reorder changes.

### Counters

However, crucially, counters are required. To this end, rather than instrument
each type with an `enum` and track it in base, we leverage CRTP to form a
registry of functions. This registry:

- Provides the calls per instance of the derived class
- Uses CRTP to ensure the count is incremented per instance as well
  - This is made clearer by the change in name from `force` to `forceImpl`

This design is closely aligned to discussions in
{cite:t}`dpt-pikusHandsonDesignPatterns2023,dpt-iglbergerSoftwareDesignDesign2022`.

We can instrument as many Potential objects as required and get both instance
counts and the number of instances.

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: DPT_
keyprefix: dpt-
---
```
