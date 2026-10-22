---
myst:
  html_meta:
    "description": "Configuration options for the Process Search job type in eOn, used within the aKMC method to find saddle points and connecting minima."
    "keywords": "eOn process search, aKMC, saddle search, minimum energy path"
---

# Process search

The aKMC method can ask clients to do a saddle search, find connecting minima,
and calculate prefactors all within the context of this job.

Those connecting minima each use a separate `ext_pot` when `[Main] parallel`
stays at its default. The exchange directories are described in
<project:ext_pot.md>.

## Known saddles

Most searches from a state end on a saddle the server already has. With
`[Process Search] known_saddles` (default `true`) the server sends each
min-mode search the saddles of its state as the job file
`known_saddles.dat`: the processes in the process table, and the saddles of
searches that converged and then failed (not connected, bad prefactor,
barrier out of window, unconverged minima, failed prefactor, negative
barrier). The failed saddles are kept in `states/<n>/knownbad/`, one
structure per distinct saddle.

At every climb point with negative curvature the client compares the
structure with each listed saddle within `[Structure Comparison]
energy_difference`, using the structure comparison. On a match the search
ends with status 23, "Climb reached a known saddle", before the two endpoint
minimizations. The server books that result as a repeat of the matching
process, or as a failed search when no process matches; a full search that
climbed to the same saddle would have ended the same way.
The key is read by both the server and the client; with `false` the server
sends no list and the client ignores one.

## Configuration

```{code-block} ini
[Process Search]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.ProcessSearchConfig
```
