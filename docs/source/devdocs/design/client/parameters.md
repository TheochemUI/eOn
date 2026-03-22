---
myst:
  html_meta:
    "description": "Design rationale for the Parameters decomposition: NSDMI defaults, INI extraction, JSON serialization, and narrow parameter passing."
    "keywords": "eOn parameters, configuration, JSON, INI, NSDMI, C++20"
---

# Parameters system

The `Parameters` class manages all runtime configuration for the eOn client.
This page documents the design decisions behind its decomposition.

## Architecture

The system has three layers:

1. **Option-group structs** (`Parameters.h`): 33 structs with C++20 NSDMI
   defaults. Each struct is independently constructable without parsing.
2. **INI parser** (`ParametersINI.cpp`): reads `config.ini` via
   [inih](https://github.com/benhoyt/inih) and populates option groups.
3. **JSON serializer** (`ParametersJSON.cpp`): provides `to_json()`/`from_json()`
   for library usage, RPC transport, and debugging.

A `validate_and_link()` function resolves cross-group dependencies (time unit
conversions, default inheritance) after loading from any source.

## NSDMI defaults

All ~290 configuration fields have default values as Non-Static Data Member
Initializers directly on the struct members:

```cpp
struct main_options_t {
  JobType job{JobType::Process_Search};
  long randomSeed{-1};
  double temperature{300.0};
  // ...
} main_options;
```

The constructor calls only `validate_and_link()` to resolve computed fields
(time unit conversions). This means a default-constructed `Parameters` object
is immediately usable.

## Cross-group dependencies

Some fields depend on values from other groups. These are resolved in
`validate_and_link()` after all groups have their raw values:

| Dependent field | Source |
|---|---|
| `optimizer_options.time_step` | `time_step_input / constants.timeUnit` |
| `neb_options.force_tolerance` | `optimizer_options.converged_force` |
| `process_search_options.minimization_offset` | `optimizer_options.max_move` |
| `dynamics_options.steps` | `floor(time / time_step)` |
| All `*_time` fields | `*_time_input / constants.timeUnit` |

## Narrow parameter passing

Core classes receive only the option groups they need via config structs:

| Class | Config struct | Fields |
|---|---|---|
| `Matter` | Direct members | `removeNetForce`, `structComp` |
| `Optimizer` hierarchy | `OptimizerConfig` | `optimizer_options_t` + 2 extras |
| `Dynamics` | `DynamicsConfig` | 12 fields from 5 groups |
| `Potential` base | `PotType` only | No Parameters dependency |

Deprecated backward-compatibility constructors are retained for callers that
still pass the full `Parameters` object.

## JSON serialization

`ParametersJSON.cpp` provides round-trip JSON serialization using
[nlohmann/json](https://github.com/nlohmann/json). This enables:

- **Library usage**: configure eOn programmatically without INI files
- **RPC transport**: send config as JSON text via capnp serve mode
- **Debugging**: dump current config to human-readable JSON

```cpp
Parameters params;
params.load_json(R"({"Main": {"job": "Nudged_Elastic_Band", "temperature": 500}})");
std::string json = params.to_json(); // pretty-printed JSON
```

## History

The v3c branch (2024) attempted to switch from INI to TOML and restructure
all parameters simultaneously. It was abandoned because it changed too many
axes at once. The current approach is incremental: extract defaults (Phase 1),
narrow passing (Phase 2), add JSON (Phase 3). Each step is independently
mergeable.
