# Guidelines

We try to follow the [C++ core
guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines).

Additionally, we fully embrace the modern value semantics of C++ over heap
oriented code.

## Preprocessor directives

- We have `EON_CHECKS` true for all but release builds
 + Meant to be used for additional checks (possibly including better errors)

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

## Maximizing Heap-Free Code

```{note}
The rationale for this is simple, `eOn`, like most scientific codes, never runs the client code for dynamic interaction.
```

We aim to minimize heap allocations and prefer stack allocation or other
techniques whenever possible.

- **Use Stack Allocation**: Prefer stack-allocated objects and automatic storage
  for local variables.
- **Use STL Containers with Stack Storage**: Utilize `std::array` for fixed-size
  arrays.
- **Avoid Dynamic Memory Allocation**: Avoid `new` and `delete`; use
  stack-allocated variables.
- **Leverage Modern C++ Features**: Use `std::optional`, `std::variant`, and
  `std::string_view` for efficient, heap-free alternatives.
