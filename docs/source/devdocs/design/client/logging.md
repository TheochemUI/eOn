# Logging

```{versionchanged} 3.x
```

`spdlog` is used for managing the logs and for writing some of the output files.
This is handled through an RAII compatible interface, since the logger is
effectively a singleton.
