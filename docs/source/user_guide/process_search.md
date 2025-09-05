---
myst:
  html_meta:
    "description": "Configuration options for the Process Search job type in EON, used within the aKMC method to find saddle points and connecting minima."
    "keywords": "EON process search, aKMC, saddle search, minimum energy path"
---

# Process Search

The aKMC method can ask clients to do a saddle search, find connecting minima,
and calculate prefactors all within the context of this job.

## Configuration

```{code-block} ini
[Process Search]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.ProcessSearchConfig
```
