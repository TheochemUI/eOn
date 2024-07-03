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

```{todo}
An `unordered_map` is read into from the TOML file, or kept as the standard TOML
output. Each class has its own parameters, which are initialized either from the
TOML output or directly through the initializer. Since we use a factory pattern
for many objects anyway, this should be fairly straightforward.
```
