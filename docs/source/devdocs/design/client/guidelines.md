# Guidelines

We try to follow the [C++ core
guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines).

## Naming

Naming is hard. Try to make it easier with these rules:

- Pure virtual classes (`ObjectiveFunction`) should have long names
- Derived classes should abbreviate, e.g. `MinObjF`

## Constructors

### Trivial constructors

For our purposes, a trivial constructor is one which only populates member
variables.

- Prefer **default member initialization**
  - Delegating constructors make DRY harder by requiring default values in
    multiple locations

## Non-trivial constructors

Try to avoid this as much as possible, factory functions should do the heavy
lifting typically, however where they are present:

- Prefer [delegating constructors](https://learn.microsoft.com/en-us/cpp/cpp/delegating-constructors?view=msvc-170)
