---
myst:
  html_meta:
    "description": "Guide to using eonclient serve mode to expose any eOn potential over RPC for integration with external tools like ChemGP."
    "keywords": "eOn serve mode, RPC server, rgpot, Cap'n Proto, ChemGP, potential serving"
---

# Serve Mode

```{versionadded} 2.2
```

```{admonition} conda-forge availability
:class: tip
Included in the `conda-forge` package (v2.12+). No additional build flags required.
```

Serve mode wraps any eOn potential as an
[rgpot](https://github.com/OmniPotentRPC/rgpot)-compatible server, exposing it
over Cap'n Proto RPC. This allows external tools, such as Julia-based
optimization frameworks (e.g.,
[ChemGP](https://github.com/HaoZeke/ChemGP)), to evaluate energies and forces
without embedding C++ code directly.

## Compilation

Serve mode requires the `with_serve` build option and a Cap'n Proto installation.
The `serve` pixi environment provides all necessary dependencies:

```{code-block} bash
pixi run -e serve bash
meson setup builddir -Dwith_serve=true
meson compile -C builddir
```

To also enable Metatomic (ML) potentials:

```{code-block} bash
meson setup builddir \
    -Dwith_serve=true \
    -Dwith_metatomic=true \
    -Dpip_metatomic=true \
    -Dtorch_version=2.9
meson compile -C builddir
```

## Usage

Serve mode supports four modes of operation: single-potential, replicated,
gateway, and multi-model. Each mode is selected by the combination of CLI flags
provided.

### Single-Potential

The simplest usage serves one potential on a single port:

```{code-block} bash
eonclient -p lj --serve-port 12345
```

This starts a blocking RPC server on `localhost:12345` serving the Lennard-Jones
potential. The server runs until interrupted with `Ctrl+C`.

To bind to all interfaces:

```{code-block} bash
eonclient -p lj --serve-host 0.0.0.0 --serve-port 12345
```

### Replicated

To start multiple copies of the same potential on sequential ports:

```{code-block} bash
eonclient -p lj --serve-port 12345 --replicas 4
```

This starts 4 independent servers on ports 12345 through 12348, each with its
own potential instance and event loop. Useful when clients can load-balance
across known ports.

### Gateway

Gateway mode exposes a single port backed by a pool of potential instances.
Incoming requests are dispatched round-robin across the pool, so clients only
need to know one address:

```{code-block} bash
eonclient -p lj --serve-port 12345 --replicas 6 --gateway
```

This creates 6 LJ potential instances and serves them all behind port 12345.
Each incoming RPC call is routed to the next available instance. This is the
recommended mode for high-throughput use cases where clients should not need to
track multiple ports.

### Multi-Model

The `--serve` flag accepts a comma-separated specification of
`potential:port` or `potential:host:port` pairs, each served concurrently
in its own thread:

```{code-block} bash
eonclient --serve "lj:12345,eam_al:12346"
```

With explicit hosts:

```{code-block} bash
eonclient --serve "lj:0.0.0.0:12345,eam_al:0.0.0.0:12346"
```

## Potential Configuration

Potentials that require parameters (Metatomic, XTB, etc.) can be configured
via the `--config` flag, which loads an INI-format config file:

```{code-block} bash
eonclient --serve "metatomic:12345" --config model.ini
```

The config file uses the same INI format as eOn's standard `config.ini`. For
example, a Metatomic model:

```{code-block} ini
[Metatomic]
model_path = /path/to/model.pt
device = cuda
length_unit = angstrom
```

Or for XTB:

```{code-block} ini
[XTBPot]
paramset = GFN2xTB
accuracy = 1.0
```

## Config-Driven Serve

The `[Serve]` section allows fully config-driven serving without CLI flags
beyond `--config`:

```{code-block} ini
[Potential]
potential = lj

[Serve]
host = localhost
port = 12345
replicas = 4
gateway_port = 0
```

```{code-block} bash
eonclient --config serve.ini
```

The dispatch logic when serving from config:

1. If `endpoints` is set, each endpoint is served in its own thread.
2. If `gateway_port > 0`, a single gateway port is opened backed by a pool
   of `replicas` potential instances.
3. Otherwise, `replicas` independent servers are started on sequential ports
   beginning at `port`.

### Examples

Gateway with a Metatomic model:

```{code-block} ini
[Potential]
potential = metatomic

[Metatomic]
model_path = /path/to/model.pt
device = cuda

[Serve]
host = 0.0.0.0
gateway_port = 12345
replicas = 6
```

Multi-model endpoints:

```{code-block} ini
[Serve]
endpoints = lj:12345,eam_al:12346,metatomic:0.0.0.0:12347

[Metatomic]
model_path = /path/to/model.pt
```

## Configuration

```{code-block} ini
[Serve]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.ServeConfig
```

## Protocol

The RPC protocol is defined by rgpot's `Potentials.capnp` schema. Each request
sends:

- **positions**: flat array `[x1, y1, z1, x2, y2, z2, ...]` (Angstroms)
- **atmnrs**: atomic numbers `[Z1, Z2, ...]`
- **box**: 3x3 cell matrix (row-major flat array)

Each response returns:

- **energy**: total potential energy (eV)
- **forces**: flat array matching positions layout (eV/Angstrom)

## Integration with ChemGP

ChemGP connects to serve mode via its `RpcPotential` oracle:

```{code-block} julia
using ChemGP

# Connect to a running eonclient --serve instance
pot = RpcPotential("localhost", 12345, atmnrs, box)
E, F = ChemGP.calculate(pot, positions)

# Use as a GP optimization oracle (gradient = -forces)
oracle = make_rpc_oracle(pot)
```

See the [ChemGP RPC tutorial](https://github.com/HaoZeke/ChemGP) for details.

## Architecture Notes

The serve mode uses a `ForceCallback` (flat-array `std::function`) interface
internally, completely decoupling the eOn potential from the capnp server.
This avoids a type collision between eOn's Eigen-based `AtomMatrix` and rgpot's
custom `AtomMatrix` by never allowing both types to coexist in the same
translation unit. The capnp schema code is compiled in a separate TU
(`ServeRpcServer.cpp`) from the eOn potential wrapper (`ServeMode.cpp`). For
more on the integration pattern, see the
[rgpot integration guide](https://rgpot.rgoswami.me/integration_guide.html).

## Command Reference

| Flag | Description |
|------|-------------|
| `--serve <spec>` | Multi-model serve spec: `pot:port` or `pot:host:port`, comma-separated |
| `--serve-host <host>` | Host for single-potential mode (default: `localhost`) |
| `--serve-port <port>` | Port for single-potential mode with `-p` (default: `12345`) |
| `--replicas <N>` | Number of server instances or gateway pool size (default: `1`) |
| `--gateway` | Enable gateway mode (single port, round-robin pool) |
| `--config <file>` | INI config file for potential and serve parameters |
| `-p <potential>` | Potential type (used with `--serve-port`) |
