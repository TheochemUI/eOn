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

## Configuration

```{code-block} ini
[Process Search]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.ProcessSearchConfig
```
