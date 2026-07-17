---
myst:
  html_meta:
    "description": "Design of eOn parameters: Cap'n Proto L0 field graph as authoring SSoT, INI/JSON adapters, and runtime Parameters."
    "keywords": "eOn parameters, Cap'n Proto, config.ini, JSON, SSoT, Parameters"
---

# Parameters system

Runtime configuration for the eOn **C++ client** lives in the `Parameters`
class. Field **names, types, and defaults** for the high-traffic option groups
are authored in a Cap'n Proto schema, not hand-synced across three peer files.

## Authoring SSoT (Cap'n Proto L0)

| Layer | Path | Role |
|-------|------|------|
| **Authoring SSoT** | {file}`schema/eon_params.capnp` | Sole place to add/rename covered fields and defaults |
| **Codegen** | `python tools/params_ssot/codegen.py` | Emits catalog + C++ defaults + field index |
| **Generated** | `schema/eon_params_catalog.json`, `eon/_params_ssot_catalog.py`, `client/generated/*` | Do not hand-edit |
| **Runtime store** | `client/Parameters.h` | Option-group structs; covered defaults applied via `apply_ssot_defaults()` |
| **Adapters** | `ParametersINI.cpp`, `ParametersJSON.cpp` | Load/save user `config.ini` / JSON |
| **Docs / validation** | `eon/schema.py`, `eon/config.yaml` | Must stay field-parity with SSoT for covered groups (`tests/test_params_ssot.py`) |

Covered groups today: **Main**, **Potential**, **Optimizer** (with LBFGS/CG/Quickmin/SD),
**Structure Comparison**, **Process Search**. Other groups still live primarily
in `Parameters.h` until folded into the schema (see `schema/README.md`).

### Adding a covered option

1. Edit {file}`schema/eon_params.capnp` (new field, type, default, ordinal).
2. Run `python tools/params_ssot/codegen.py`.
3. Wire any INI/JSON load path if needed; do **not** invent the same field only
   in `config.yaml` / `schema.py` / `Parameters.h` without the Cap'n Proto update.

Users still write ordinary **`config.ini`** (or JSON) files. Cap'n Proto is the
**field-graph author**, not a requirement that end users ship binary messages.

## Runtime architecture

1. **Option-group structs** (`Parameters.h`): nested C++ structs; covered
   defaults originate from the Cap'n Proto schema via
   `eonc::config::apply_ssot_defaults()` in the `Parameters` constructor.
2. **INI parser** (`ParametersINI.cpp`): reads `config.ini` via
   [inih](https://github.com/benhoyt/inih) and populates option groups.
3. **JSON serializer** (`ParametersJSON.cpp`): `to_json()` / `from_json()` for
   library use, RPC text blobs, and debugging.

`validate_and_link()` resolves cross-group dependencies (time unit conversions,
default inheritance) after loading from any source.

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

## Python server vs client

The **Python server** (`python -m eon` / `eon-server`) still uses
`eon.config.ConfigClass` driven by `config.ini` and `eon/config.yaml` for
validation of allowed keys. For covered groups, those keys and defaults are
parity-checked against the Cap'n Proto catalog. The client binary
(`eonclient`) uses `Parameters` with the same SSoT-backed defaults for those
groups.

## History

The v3c branch (2024) attempted to switch from INI to TOML and restructure
all parameters simultaneously. It was abandoned because it changed too many
axes at once. Later work added NSDMI defaults, JSON I/O, and narrow passing.
The Cap'n Proto L0 field graph unifies **authoring** for the core groups without
forcing a pure-binary config format on users.
