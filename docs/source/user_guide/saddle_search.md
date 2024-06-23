# Saddle Search

A saddle search is initiated by making a local displacement of atoms from their
position at the minimum of the current state. This displacement can be done
using the different strategies indicated by the {any}`client_displace_type`
option, and the following parameters. If the user knows something about the
local environment where reactions are likely to take place in the system, this
information can be used to make saddle searches more efficient by getting them
started in the right part of configuration space.

## Configuration

```{code-block} ini
[Saddle Search]
```


```{versionchanged} 2.1_TBA
In TOML, this will be `[Saddle_Search]`
```


```{eval-rst}
.. autopydantic_model:: eon.schema.SaddleSearchConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: SS_
keyprefix: ss-
---
```
