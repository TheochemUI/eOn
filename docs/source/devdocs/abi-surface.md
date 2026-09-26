# Public C++ surface

The installed `include/eon` tree still ships internal headers. The
Job factory has landed. `Parameters` stores load state and option groups
in a private `Impl`, so `sizeof(Parameters)` is that pointer. Const
accessors and `ParametersLoadAccess` are the read and write surface.
Option-group types remain in the installed `ParametersOptions.h`.
`Matter` still exposes Eigen members in the header (accessors return
Eigen types). Treat only these as the extension contract:

| Type | Pure virtuals | Notes |
|---|---|---|
| `Potential` | `force(...)` | Raw C arrays stay the Fortran/FFI adapter. `get_ef` is Eigen. |
| `Job` | `run()` | `runFromMatter` is extra, not required. |
| `Optimizer` | `step`, `run` | New methods must be defaulted virtuals, not new pures. |

Capability queries (`isThreadSafe`, `layoutFlags`, `supportsBatchEvaluation`)
are defaulted virtuals. Do not add another `= 0` without a towncrier
breaking fragment.

CI should keep grepping that those three pures still exist and that new
pures are not added silently.
